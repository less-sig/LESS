// #[inline(always)]
// pub fn compute_ct_mask(i: usize, j: usize) -> isize {
//     -((i == j) as isize)
// }
use num_traits::{
    PrimInt, WrappingNeg
};

#[inline(always)]
pub fn compute_ct_mask<T>(i: T, j: T) -> T
where
    T: PrimInt + WrappingNeg,
{
    let t: bool = i == j;
    T::from(t as u8).unwrap().wrapping_neg()
}
