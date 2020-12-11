mod cubic_poly;
mod zero;

use std::{iter, ops};

pub use cubic_poly::CubicPoly;
pub use zero::Zero;

/// Represents boundary conditions to be used for fitting a spline
pub enum BoundaryCondition<T> {
    /// Set derivatives at the initial and final points
    Derivatives(T, T),
    /// Set second derivatives at the initial and final points
    SecondDerivatives(T, T),
    /// Second derivatives at initial and final points set to 0
    /// (equivalent to SecondDerivatives(0.0, 0.0))
    Natural,
    /// fit a periodic function
    Periodic,
}

/// A result of interpolation between a set of points
#[derive(Clone, Debug)]
pub struct Spline<T> {
    points_x: Vec<f64>,
    splines: Vec<CubicPoly<T>>,
    derivative_start: T,
    derivative_end: T,
}

impl<T> Spline<T>
where
    T: ops::Add<T, Output = T>
        + ops::AddAssign<T>
        + ops::Sub<T, Output = T>
        + ops::SubAssign<T>
        + ops::Mul<f64, Output = T>
        + ops::Div<f64, Output = T>
        + Copy
        + Zero,
{
    fn with_derivatives(points: Vec<(f64, T)>, d0: T, dn: T) -> Self {
        let n = points.len();
        let points_x: Vec<f64> = points.iter().map(|&(x, _)| x).collect();
        let points_y: Vec<T> = points.iter().map(|&(_, y)| y).collect();
        let h: Vec<f64> = points_x.windows(2).map(|w| w[1] - w[0]).collect();

        let a: Vec<f64> = iter::once(0.0)
            .chain(h.windows(2).map(|w| w[0] / (w[0] + w[1])))
            .chain(iter::once(1.0))
            .collect();
        let b: Vec<f64> = iter::repeat(2.0).take(points_x.len()).collect();
        let c: Vec<f64> = iter::once(1.0)
            .chain(h.windows(2).map(|w| w[1] / (w[0] + w[1])))
            .chain(iter::once(0.0))
            .collect();
        let d: Vec<T> = iter::once((div_diff_2(points[0], points[1]) - d0) * 6.0 / h[0])
            .chain(
                points
                    .windows(3)
                    .map(|w| div_diff_3(w[0], w[1], w[2]) * 6.0),
            )
            .chain(iter::once(
                (dn - div_diff_2(points[n - 2], points[n - 1])) * 6.0 / h[n - 2],
            ))
            .collect();
        let second_derivatives = solve_tridiagonal(a, b, c, d);
        let splines = (0..n - 1)
            .map(|i| {
                let hi = h[i];
                let mi = second_derivatives[i];
                let mi1 = second_derivatives[i + 1];
                let yi = points_y[i];
                let yi1 = points_y[i + 1];
                let poly1 = CubicPoly::new(
                    mi / -6.0 / hi,
                    T::zero(),
                    (yi - mi * hi * hi / 6.0) / -hi,
                    T::zero(),
                )
                .shifted(points_x[i + 1]);
                let poly2 = CubicPoly::new(
                    mi1 / 6.0 / hi,
                    T::zero(),
                    (yi1 - mi1 * hi * hi / 6.0) / hi,
                    T::zero(),
                )
                .shifted(points_x[i]);
                poly1 + poly2
            })
            .collect();
        Self {
            points_x,
            splines,
            derivative_start: d0,
            derivative_end: dn,
        }
    }

    fn with_second_derivatives(points: Vec<(f64, T)>, d0: T, dn: T) -> Self {
        let n = points.len();
        let points_x: Vec<f64> = points.iter().map(|&(x, _)| x).collect();
        let points_y: Vec<T> = points.iter().map(|&(_, y)| y).collect();
        let h: Vec<f64> = points_x.windows(2).map(|w| w[1] - w[0]).collect();

        let a: Vec<f64> = iter::once(0.0)
            .chain(h.windows(2).map(|w| w[0] / (w[0] + w[1])))
            .chain(iter::once(0.0))
            .collect();
        let b: Vec<f64> = iter::repeat(2.0).take(points_x.len()).collect();
        let c: Vec<f64> = iter::once(0.0)
            .chain(h.windows(2).map(|w| w[1] / (w[0] + w[1])))
            .chain(iter::once(0.0))
            .collect();
        let d: Vec<T> = iter::once(d0 * 2.0)
            .chain(
                points
                    .windows(3)
                    .map(|w| div_diff_3(w[0], w[1], w[2]) * 6.0),
            )
            .chain(iter::once(dn * 2.0))
            .collect();
        let second_derivatives = solve_tridiagonal(a, b, c, d);
        let splines: Vec<_> = (0..n - 1)
            .map(|i| {
                let hi = h[i];
                let mi = second_derivatives[i];
                let mi1 = second_derivatives[i + 1];
                let yi = points_y[i];
                let yi1 = points_y[i + 1];
                let poly1 = CubicPoly::new(
                    mi / -6.0 / hi,
                    T::zero(),
                    (yi - mi * hi * hi / 6.0) / -hi,
                    T::zero(),
                )
                .shifted(points_x[i + 1]);
                let poly2 = CubicPoly::new(
                    mi1 / 6.0 / hi,
                    T::zero(),
                    (yi1 - mi1 * hi * hi / 6.0) / hi,
                    T::zero(),
                )
                .shifted(points_x[i]);
                poly1 + poly2
            })
            .collect();
        let x0 = points_x[0];
        let xn = points_x[n - 1];
        let d0 = splines[0].derivative(x0);
        let dn = splines[n - 2].derivative(xn);
        Self {
            points_x,
            splines,
            derivative_start: d0,
            derivative_end: dn,
        }
    }

    /// Creates a new interpolated function fit to a set of given points with the given boundary
    /// conditions
    pub fn new(mut points: Vec<(f64, T)>, boundary_condition: BoundaryCondition<T>) -> Self {
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        match boundary_condition {
            BoundaryCondition::Derivatives(d0, dn) => Self::with_derivatives(points, d0, dn),
            BoundaryCondition::SecondDerivatives(d0, dn) => {
                Self::with_second_derivatives(points, d0, dn)
            }
            BoundaryCondition::Natural => {
                Self::with_second_derivatives(points, T::zero(), T::zero())
            }
            BoundaryCondition::Periodic => unimplemented!(),
        }
    }

    /// Evaluates the interpolated function at a given point
    pub fn eval(&self, x: f64) -> T {
        let index = self
            .points_x
            .binary_search_by(|probe| probe.partial_cmp(&x).unwrap());
        match index {
            Ok(i) => {
                if i < self.splines.len() {
                    self.splines[i].eval(x)
                } else {
                    self.splines[i - 1].eval(x)
                }
            }
            Err(i) => {
                if i == 0 {
                    let x0 = self.points_x[0];
                    let y0 = self.splines[0].eval(x0);
                    y0 + self.derivative_start * (x - x0)
                } else if i == self.points_x.len() {
                    let xn = self.points_x[self.points_x.len() - 1];
                    let yn = self.splines[self.splines.len() - 1].eval(xn);
                    yn + self.derivative_end * (x - xn)
                } else {
                    self.splines[i - 1].eval(x)
                }
            }
        }
    }

    pub fn eval_derivative(&self, x: f64) -> T {
        let index = self
            .points_x
            .binary_search_by(|probe| probe.partial_cmp(&x).unwrap());
        match index {
            Ok(i) => {
                if i < self.splines.len() {
                    self.splines[i].derivative(x)
                } else {
                    self.splines[i - 1].derivative(x)
                }
            }
            Err(i) => {
                if i == 0 {
                    self.derivative_start
                } else if i == self.points_x.len() {
                    self.derivative_end
                } else {
                    self.splines[i - 1].derivative(x)
                }
            }
        }
    }

    pub fn derivative_start(&self) -> T {
        self.derivative_start
    }

    pub fn derivative_end(&self) -> T {
        self.derivative_end
    }

    pub fn polynomials<'a>(&'a self) -> impl Iterator<Item = (f64, f64, CubicPoly<T>)> + 'a {
        self.points_x
            .windows(2)
            .zip(self.splines.iter())
            .map(|(xs, poly)| (xs[0], xs[1], *poly))
    }
}

fn div_diff_2<T>((x0, y0): (f64, T), (x1, y1): (f64, T)) -> T
where
    T: ops::Sub<T, Output = T> + ops::Div<f64, Output = T>,
{
    (y1 - y0) / (x1 - x0)
}

fn div_diff_3<T>((x0, y0): (f64, T), (x1, y1): (f64, T), (x2, y2): (f64, T)) -> T
where
    T: ops::Add<T, Output = T> + ops::Div<f64, Output = T>,
{
    y0 / (x0 - x1) / (x0 - x2) + y1 / (x1 - x0) / (x1 - x2) + y2 / (x2 - x0) / (x2 - x1)
}

fn solve_tridiagonal<T>(a: Vec<f64>, mut b: Vec<f64>, c: Vec<f64>, mut d: Vec<T>) -> Vec<T>
where
    T: ops::Sub<T, Output = T>
        + ops::SubAssign<T>
        + ops::Mul<f64, Output = T>
        + ops::Div<f64, Output = T>
        + Copy,
{
    let n = b.len();
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= c[i - 1] * w;
        let z = d[i - 1] * w;
        d[i] -= z;
    }
    let mut result = vec![d[n - 1] / b[n - 1]];
    for i in (0..n - 1).rev() {
        let last_x = result[result.len() - 1];
        result.push((d[i] - last_x * c[i]) / b[i]);
    }
    result.reverse();
    result
}

#[cfg(test)]
mod tests {
    use super::{BoundaryCondition, CubicPoly, Spline};

    #[test]
    fn test_spline_with_derivatives() {
        let points = vec![
            (0.0, 0.0),
            (1.0, 6.0),
            (1.2, 6.0),
            (1.4, 6.0),
            (1.6, 6.0),
            (2.0, 1.0),
            (3.0, 2.0),
            (4.0, -1.0),
        ];
        let spline = Spline::new(points, BoundaryCondition::Natural);
        for i in -20..220 {
            let x = (i as f64) / 200.0 * 4.0;
            println!("{} {}", x, spline.eval(x));
        }
    }

    #[test]
    fn test_atmosphere() {
        let points = vec![
            (0.0, 8.6),
            (12.5, 10.4),
            (19.4, 8.9),
            (24.0, 11.7),
            (34.0, 17.5),
        ];
        let spline = Spline::new(points, BoundaryCondition::Derivatives(-0.0065, -0.0065));
        println!("{:?}", spline);
        for i in -20..240 {
            let x = (i as f64) * 0.5;
            println!("{} {}", x, spline.eval(x));
        }
    }

    #[test]
    fn test_spline_of_polynomials() {
        let points = vec![
            (0.0, CubicPoly::new(1.0, 0.0, -1.0, 0.0)),
            (1.0, CubicPoly::new(2.0, -0.5, -3.0, 2.4)),
            (2.0, CubicPoly::new(3.0, -1.5, 1.0, -1.2)),
            (3.0, CubicPoly::new(2.0, 1.5, 2.0, 1.5)),
            (4.0, CubicPoly::new(1.0, 3.0, -3.0, -2.0)),
        ];
        let spline = Spline::new(points, BoundaryCondition::Natural);
        for i in -10..110 {
            let x = (i as f64) / 100.0 * 4.0;
            let y = 8.0;
            println!("{} {} {}", x, y, spline.eval(x).eval(y));
        }
    }
}
