pub struct CubicPoly {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

impl CubicPoly {
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }

    /// Creates a polynomial g(x) = f(x-x0), where f(x) = self
    pub fn shifted(self, x0: f64) -> Self {
        let a = self.a;
        let b = self.b - 3.0 * self.a * x0;
        let c = self.c + 3.0 * self.a * x0 * x0 - 2.0 * self.b * x0;
        let d = self.d - self.a * x0 * x0 * x0 + self.b * x0 * x0 - self.c * x0;
        Self { a, b, c, d }
    }

    pub fn eval(&self, x: f64) -> f64 {
        self.a * x * x * x + self.b * x * x + self.c * x + self.d
    }
}

pub enum BoundaryCondition {
    Derivatives(f64, f64),
    SecondDerivatives(f64, f64),
    Natural,
    Periodic,
}

pub struct Spline {
    points_x: Vec<f64>,
    splines: Vec<CubicPoly>,
    derivative_start: f64,
    derivative_end: f64,
}

impl Spline {
    pub fn new(points: Vec<(f64, f64)>, boundary_condition: BoundaryCondition) -> Self {
        // TODO
        Self {
            points_x: Default::default(),
            splines: Default::default(),
            derivative_start: 0.0,
            derivative_end: 0.0,
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
