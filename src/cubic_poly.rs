use std::ops;

#[derive(Clone, Copy, Debug)]
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

    pub fn derivative(&self, x: f64) -> f64 {
        3.0 * self.a * x * x + 2.0 * self.b * x + self.c
    }
}

impl ops::AddAssign<CubicPoly> for CubicPoly {
    fn add_assign(&mut self, other: CubicPoly) {
        self.a += other.a;
        self.b += other.b;
        self.c += other.c;
        self.d += other.d;
    }
}

impl ops::SubAssign<CubicPoly> for CubicPoly {
    fn sub_assign(&mut self, other: CubicPoly) {
        self.a -= other.a;
        self.b -= other.b;
        self.c -= other.c;
        self.d -= other.d;
    }
}

impl ops::Add<CubicPoly> for CubicPoly {
    type Output = CubicPoly;

    fn add(mut self, other: CubicPoly) -> CubicPoly {
        self += other;
        self
    }
}

impl ops::Sub<CubicPoly> for CubicPoly {
    type Output = CubicPoly;

    fn sub(mut self, other: CubicPoly) -> CubicPoly {
        self -= other;
        self
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
