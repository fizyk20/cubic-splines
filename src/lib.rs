mod cubic_poly;

pub use cubic_poly::CubicPoly;
use std::iter;

pub enum BoundaryCondition {
    Derivatives(f64, f64),
    SecondDerivatives(f64, f64),
    Natural,
    Periodic,
}

#[derive(Clone, Debug)]
pub struct Spline {
    points_x: Vec<f64>,
    splines: Vec<CubicPoly>,
    derivative_start: f64,
    derivative_end: f64,
}

impl Spline {
    fn with_derivatives(points: Vec<(f64, f64)>, d0: f64, dn: f64) -> Self {
        let n = points.len();
        let points_x: Vec<f64> = points.iter().map(|&(x, _)| x).collect();
        let points_y: Vec<f64> = points.iter().map(|&(_, y)| y).collect();
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
        let d: Vec<f64> = iter::once(6.0 / h[0] * (div_diff_2(points[0], points[1]) - d0))
            .chain(
                points
                    .windows(3)
                    .map(|w| 6.0 * div_diff_3(w[0], w[1], w[2])),
            )
            .chain(iter::once(
                6.0 / h[n - 2] * (dn - div_diff_2(points[n - 2], points[n - 1])),
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
                let poly1 =
                    CubicPoly::new(-mi / 6.0 / hi, 0.0, -(yi - mi * hi * hi / 6.0) / hi, 0.0)
                        .shifted(points_x[i + 1]);
                let poly2 =
                    CubicPoly::new(mi1 / 6.0 / hi, 0.0, (yi1 - mi1 * hi * hi / 6.0) / hi, 0.0)
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

    pub fn new(mut points: Vec<(f64, f64)>, boundary_condition: BoundaryCondition) -> Self {
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        match boundary_condition {
            BoundaryCondition::Derivatives(d0, dn) => Self::with_derivatives(points, d0, dn),
            BoundaryCondition::SecondDerivatives(_d0, _dn) => unimplemented!(),
            BoundaryCondition::Natural => unimplemented!(),
            BoundaryCondition::Periodic => unimplemented!(),
        }
    }

    pub fn eval(&self, x: f64) -> f64 {
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
}

fn div_diff_2((x0, y0): (f64, f64), (x1, y1): (f64, f64)) -> f64 {
    (y1 - y0) / (x1 - x0)
}

fn div_diff_3((x0, y0): (f64, f64), (x1, y1): (f64, f64), (x2, y2): (f64, f64)) -> f64 {
    y0 / (x0 - x1) / (x0 - x2) + y1 / (x1 - x0) / (x1 - x2) + y2 / (x2 - x0) / (x2 - x1)
}

fn solve_tridiagonal(a: Vec<f64>, mut b: Vec<f64>, c: Vec<f64>, mut d: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut result = vec![d[n - 1] / b[n - 1]];
    for i in (0..n - 1).rev() {
        let last_x = result[result.len() - 1];
        result.push((d[i] - c[i] * last_x) / b[i]);
    }
    result.reverse();
    result
}

#[cfg(test)]
mod tests {
    use super::CubicPoly;

    #[test]
    fn test_poly_shift() {
        let poly = CubicPoly::new(1.0, -1.0, 1.0, -1.0);
        assert_eq!(poly.eval(0.0), -1.0);
        assert_eq!(poly.eval(1.0), 0.0);
        assert_eq!(poly.eval(2.0), 5.0);
        let poly2 = poly.shifted(1.0); // poly2(x) = poly(x - 1)
        assert_eq!(poly2.eval(1.0), -1.0);
        assert_eq!(poly2.eval(2.0), 0.0);
        assert_eq!(poly2.eval(3.0), 5.0);
    }
}
