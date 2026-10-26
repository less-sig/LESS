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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_ct_mask_u8() {
        for i in 0..=255u8 {
            for j in 0..=255u8 {
                let m = compute_ct_mask(i, j);
                assert_eq!(m, if i == j { 0xffu8 } else { 0x00u8 });
            }
        }
    }

    #[test]
    fn compute_ct_mask_u16() {
        assert_eq!(compute_ct_mask(0u16, 0u16), 0xffffu16);
        assert_eq!(compute_ct_mask(1234u16, 1234u16), 0xffffu16);
        assert_eq!(compute_ct_mask(1234u16, 1235u16), 0u16);
        assert_eq!(compute_ct_mask(0u16, 0xffffu16), 0u16);
    }

    #[test]
    fn compute_ct_mask_u32() {
        assert_eq!(compute_ct_mask(0u32, 0u32), 0xffff_ffffu32);
        assert_eq!(compute_ct_mask(7u32, 8u32), 0u32);
    }

    /// the mask is used to select between two values without branching
    #[test]
    fn compute_ct_mask_select() {
        for i in 0..64u16 {
            let mask = compute_ct_mask(i, 17u16);
            let v = (0xaaaau16 & !mask) ^ (0x5555u16 & mask);
            assert_eq!(v, if i == 17 { 0x5555 } else { 0xaaaa });
        }
    }
}
