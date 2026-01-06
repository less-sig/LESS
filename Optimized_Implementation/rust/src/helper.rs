// #[inline(always)]
// pub fn compute_ct_mask(i: usize, j: usize) -> isize {
//     -((i == j) as isize)
// }

#[inline(always)]
pub fn compute_ct_mask<T>(i: T, j: T) -> T
where
    T: PartialEq,
{
    T::from(i == j).unwrap().wrapping_neg()
}
