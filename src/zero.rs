use std::ops;

use super::CubicPoly;

pub trait Zero {
    fn zero() -> Self;
}

impl Zero for f64 {
    fn zero() -> f64 {
        0.0
    }
}

impl<T> Zero for CubicPoly<T>
where
    T: ops::Add<T, Output = T>
        + ops::Sub<T, Output = T>
        + ops::AddAssign<T>
        + ops::SubAssign<T>
        + ops::Mul<f64, Output = T>
        + Zero
        + Copy,
{
    fn zero() -> CubicPoly<T> {
        CubicPoly::new(T::zero(), T::zero(), T::zero(), T::zero())
    }
}
