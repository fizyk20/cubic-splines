use roots::{find_roots_cubic, Roots};
use std::ops;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicPoly {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Factors {
    /// f(x) = a(x-x1)(x-x2)(x-x3)
    ThreeLinear { a: f64, x1: f64, x2: f64, x3: f64 },
    /// f(x) = a(x-x1)(x²+bx+c)
    LinearAndQuadratic { a: f64, x1: f64, b: f64, c: f64 },
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

    pub fn factors(&self) -> Factors {
        let roots = find_roots_cubic(self.a, self.b, self.c, self.d);
        match roots {
            Roots::One([x1]) | Roots::Two([x1, _]) => {
                let b = self.b / self.a + x1;
                let c = self.c / self.a + b * x1;
                // make sure that we haven't missed any real roots
                let delta = b * b - 4.0 * c;
                if delta >= 0.0 {
                    let x2 = 0.5 * (-b - delta.sqrt());
                    let x3 = 0.5 * (-b + delta.sqrt());
                    // sort the roots
                    let (x1, x2) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
                    let (x1, x3) = if x1 < x3 { (x1, x3) } else { (x3, x1) };
                    let (x2, x3) = if x2 < x3 { (x2, x3) } else { (x3, x2) };
                    Factors::ThreeLinear {
                        a: self.a,
                        x1,
                        x2,
                        x3,
                    }
                } else {
                    Factors::LinearAndQuadratic {
                        a: self.a,
                        x1,
                        b,
                        c,
                    }
                }
            }
            Roots::Three([x1, x2, x3]) => Factors::ThreeLinear {
                a: self.a,
                x1,
                x2,
                x3,
            },
            _ => panic!("should have either one or three roots! {:?}", roots),
        }
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
    use super::{CubicPoly, Factors};

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

    #[test]
    fn test_triple_root() {
        let poly = CubicPoly::new(2.0, -6.0, 6.0, -2.0);
        assert_eq!(
            poly.factors(),
            Factors::ThreeLinear {
                a: 2.0,
                x1: 1.0,
                x2: 1.0,
                x3: 1.0,
            }
        );
    }

    #[test]
    fn test_double_root() {
        let poly = CubicPoly::new(1.0, 1.0, -1.0, -1.0);
        assert_eq!(
            poly.factors(),
            Factors::ThreeLinear {
                a: 1.0,
                x1: -1.0,
                x2: -1.0,
                x3: 1.0,
            }
        );
    }

    #[test]
    fn test_single_root() {
        let poly = CubicPoly::new(1.0, -1.0, 1.0, -1.0);
        assert_eq!(
            poly.factors(),
            Factors::LinearAndQuadratic {
                a: 1.0,
                x1: 1.0,
                b: 0.0,
                c: 1.0,
            }
        );
    }
}
