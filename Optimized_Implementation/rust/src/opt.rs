#![allow(unused_imports)]

use std::arch::x86_64::{
    __m256i,
    _mm256_set1_epi8, _mm256_extract_epi32,
    _mm256_add_epi8, _mm256_sub_epi8, _mm256_min_epu8, _mm256_blendv_epi8
};

/// using blend instruction
/// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let b = _mm256_set1_epi8(j as i8);
/// let c = gf127_add_u256(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
fn gf127_add_u256(a: __m256i, b: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_add_epi8(a, b);
    let t2 = _mm256_sub_epi8(t1, q);
    let t3 = _mm256_blendv_epi8(t2, t1, t2);
    t3
}

/// using the min instruction
/// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let b = _mm256_set1_epi8(j as i8);
/// let c = gf127_add_u256(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
fn gf127_add_u256_v2(a: __m256i, b: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_add_epi8(a, b);
    let t2 = _mm256_sub_epi8(t1, q);
    let t3 = _mm256_min_epu8(t1, t2);
    t3
}

/// using the min instruction
/// c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let b = _mm256_set1_epi8(j as i8);
/// let c = gf127_sub_u256(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_sub_u256(a: __m256i, b: __m256i) -> __m256i {
    let q  = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_sub_epi8(a, b);
    let t2 = _mm256_add_epi8(t1, q);
    let t3 = _mm256_blendv_epi8(t1, t2, t1);
    t3
} 

/// using the min instruction
/// c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let b = _mm256_set1_epi8(j as i8);
/// let c = gf127_sub_u256(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_sub_u256_v2(a: __m256i, b: __m256i) -> __m256i {
    let q  = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_sub_epi8(a, b);
    let t2 = _mm256_add_epi8(t1, q);
    let t3 = _mm256_min_epu8(t1, t2);
    t3
} 


// #[cfg(target_feature = "avx2")]
// unsafe fn gf127_mul_u256(a1: __m256i, a2: __m256i,
//                          b1: __m256i, b2: __m256i) -> __m256i {
//     let c7f = _mm256_set1_epi8(127);
//     let c007f = _mm256_set1_epi16(127);
// 
//     let mut a_lo: __m256i = a1;
//     let mut a_hi: __m256i = a2;
//     let b_lo: __m256i = b1;
//     let b_hi: __m256i = b2;
// 
//     a_lo = _mm256_mullo_epi16(a_lo, b_lo);
//     a_hi = _mm256_mullo_epi16(a_hi, b_hi);
// 
//     let mut v1 = _mm256_srli_epi16(a_lo, 7);
//     let mut v2 = _mm256_srli_epi16(a_hi, 7);
// 
//     let mut w1 = _mm256_packs_epi16(_mm256_and_si256(a_lo, c007f), _mm256_and_si256(a_hi, c007f));
//     let mut w2 = _mm256_packs_epi16(v1, v2);
// 
//     v1 = _mm256_add_epi8(w1, w2);
//     v2 = _mm256_permute4x64_epi64(v1, 0b11011000);
//     w1 = _mm256_sub_epi8(v2, c7f);
//     // w2 = _mm256_blendv_epi8(w1, v2, w1);
//     w2 = _mm256_min_epu8(v2, w1);
//     return w2;
// }


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf127_add_u256() { 
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm256_set1_epi8(i as i8);
                    let b = _mm256_set1_epi8(j as i8);
                    let c = gf127_add_u256(a, b);

                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, (i+j) % 127);
                }
            }
        }
    }

    #[test]
    fn test_gf127_add_u256_v2() { 
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm256_set1_epi8(i as i8);
                    let b = _mm256_set1_epi8(j as i8);
                    let c = gf127_add_u256_v2(a, b);

                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, (i+j) % 127);
                }
            }
        }
    }

    #[test]
    fn test_gf127_sub_u256() { 
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm256_set1_epi8(i as i8);
                    let b = _mm256_set1_epi8(j as i8);
                    let c = gf127_sub_u256(a, b);

                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, (i+127-j) % 127);
                }
            }
        }
    }

    #[test]
    fn test_gf127_sub_u256_v2() { 
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm256_set1_epi8(i as i8);
                    let b = _mm256_set1_epi8(j as i8);
                    let c = gf127_sub_u256_v2(a, b);

                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, (i+127-j) % 127);
                }
            }
        }
    }
}
