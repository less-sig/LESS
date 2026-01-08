#![allow(unused_imports)]

use crate::{fq::Fq, vector::Vector};

use std::arch::x86_64::{
    __m128i, __m256i,
    _mm_loadu_si128, _mm256_loadu_si256,
    _mm_storeu_si128, _mm256_storeu_si256,
    _mm256_setzero_si256, _mm_set1_epi8,_mm256_set1_epi8, _mm256_set1_epi16, _mm256_extract_epi32, _mm256_extract_epi8, 
    _mm256_add_epi8, _mm256_sub_epi8, _mm256_min_epu8, 
    _mm256_permute4x64_epi64, _mm256_permute2x128_si256, _mm256_blendv_epi8,
    _mm256_srli_epi16, _mm256_srli_epi32, _mm256_srli_epi64, _mm256_srli_si256,
    _mm256_mullo_epi16, _mm256_cvtepu8_epi16, 
    _mm256_packs_epi16, _mm256_and_si256
};

/// using blend instruction
/// c[i] = a[i] mod q,  a[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let c = gf127_red_u256(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
#[target_feature(enable = "avx,avx2")]
fn gf127_red_u256(a: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let b = _mm256_sub_epi8(a, q);
    let t3 = _mm256_blendv_epi8(b, a, b);
    t3
}

/// using min instruction
/// c[i] = a[i] mod q,  a[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let c = gf127_red_u256(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
#[target_feature(enable = "avx,avx2")]
fn gf127_red_u256_v2(a: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let b = _mm256_sub_epi8(a, q);
    let t3 = _mm256_min_epu8(b, a);
    t3
}

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
/// let a = _mm256_set1_epi8(0 as i8);
/// let b = _mm256_set1_epi8(1 as i8);
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


/// using the min instruction
///  sum_{i=0}^{i < 32} v[i] mod q, v[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(1 as i8);
/// let c = gf127_hadd_avx2(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_hadd_avx2(v: __m256i) -> u8 {
    let mut a: __m256i  = _mm256_srli_epi16(v, 8);
    let mut t: __m256i  = _mm256_add_epi8(a, v);
    t = gf127_red_u256(t);

    a = _mm256_srli_epi32(t, 16);
    t = _mm256_add_epi8(a, t);
    t = gf127_red_u256(t);

    a = _mm256_srli_epi64(t, 32);
    t = _mm256_add_epi8(a, t);
    t = gf127_red_u256(t);

    a = _mm256_srli_si256(t, 8);
    t = _mm256_add_epi8(a, t);
    t = gf127_red_u256(t);

    a = _mm256_permute2x128_si256(t, t, 1);
    t = _mm256_add_epi8(a, t);
    t = gf127_red_u256(t);

    return _mm256_extract_epi8(t, 0) as u8;
}

/// TODO test
/// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(0 as i8);
/// let b = _mm256_set1_epi8(1 as i8);
/// let c = gf127_sub_u256(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_scalar_mul_avx2(a: [Fq; 32], b: __m256i) -> __m256i {
    let c7f = _mm256_set1_epi8(127);
    let c007f = _mm256_set1_epi16(127);

    let t1 = _mm_loadu_si128(a.as_ptr() as *const __m128i);
    let t2 = _mm_loadu_si128(a.as_ptr().add(16) as *const __m128i);
    
    let mut a_lo: __m256i = _mm256_cvtepu8_epi16(t1);
    let mut a_hi: __m256i = _mm256_cvtepu8_epi16(t2);
    
    a_lo = _mm256_mullo_epi16(a_lo, b);
    a_hi = _mm256_mullo_epi16(a_hi, b);
    
    let mut v1 = _mm256_srli_epi16(a_lo, 7);
    let mut v2 = _mm256_srli_epi16(a_hi, 7);

    let w1: __m256i = _mm256_packs_epi16(_mm256_and_si256(a_lo, c007f), _mm256_and_si256(a_hi, c007f));
    let w2: __m256i = _mm256_packs_epi16(v1, v2);

    v1 = _mm256_add_epi8(w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0xd8); // 0b11011000
    let t = gf127_red_u256(v2);
    return t;
}

/// c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a1 = _mm_set1_epi8(0 as i8);
/// let a2 = _mm_set1_epi8(0 as i8);
/// let b1 = _mm_set1_epi8(0 as i8);
/// let b2 = _mm_set1_epi8(0 as i8);
/// let c = gf127_mul_avx2(a1, a2, b1, b2);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_mul_avx2(a1: __m128i, a2: __m128i,
                         b1: __m128i, b2: __m128i) -> __m256i {
    let c007f = _mm256_set1_epi16(127);

    let mut a_lo: __m256i = _mm256_cvtepu8_epi16(a1);
    let mut a_hi: __m256i = _mm256_cvtepu8_epi16(a2);
    let b_lo: __m256i = _mm256_cvtepu8_epi16(b1);
    let b_hi: __m256i = _mm256_cvtepu8_epi16(b2);

    a_lo = _mm256_mullo_epi16(a_lo, b_lo);
    a_hi = _mm256_mullo_epi16(a_hi, b_hi);
    
    let mut v1 = _mm256_srli_epi16(a_lo, 7);
    let mut v2 = _mm256_srli_epi16(a_hi, 7);

    let w1: __m256i = _mm256_packs_epi16(_mm256_and_si256(a_lo, c007f), _mm256_and_si256(a_hi, c007f));
    let w2: __m256i = _mm256_packs_epi16(v1, v2);

    v1 = _mm256_add_epi8(w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0xd8); // 0b11011000
    let t = gf127_red_u256(v2);
    return t;
}

/// \return sum_i..N a[i] mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let a1 = gf127_row_acc_avx2(&mut row);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_acc_avx2<const N: usize>(row: &mut Vector<N>) -> u8 {
    assert!(N % 32 == 0);
    unsafe {
        let ptr = row.as_ptr(); // *const u8
        let mut s: __m256i = _mm256_loadu_si256(ptr as *const __m256i);

        for col in (32..N).step_by(32) {
            let t1: __m256i = _mm256_loadu_si256(ptr.add(col) as *const __m256i);
            s = gf127_add_u256(s, t1); 
        }

        return gf127_hadd_avx2(s);
    }
}

/// \return sum_i..N a[i]^-1 mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let a1 = gf127_row_acc_inv_avx2(&mut row);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_acc_inv_avx2<const N: usize>(row: &mut Vector<N>) -> u8 {
    assert!(N % 32 == 0);
    let mut inv_data: Vector<N> = Vector::new();
    for col in 0..N {
        inv_data[col] = Fq::inv_non_ct(row[col]);
    }

    return gf127_row_acc_avx2(&mut inv_data);
}

/// \return a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_avx2(&mut row), s;
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_scalar_mul_avx2<const N: usize>(row: &mut Vector<N>,
                                                 s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr = row.as_ptr(); // *const u8
        let b = _mm_set1_epi8(s.0 as i8);

        for col in (0..N).step_by(32) {
            let a1: __m128i = _mm_loadu_si128(ptr.add(col +  0) as *const __m128i);
            let a2: __m128i = _mm_loadu_si128(ptr.add(col + 16) as *const __m128i);
            let t2: __m256i = gf127_mul_avx2(a1, a2, b, b); 

            _mm256_storeu_si256(ptr.add(col) as *mut __m256i, t2);
        }
    }
}

/// c[i] = a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row_out = Vector::<N>::new();
/// let row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_2_avx2(&mut row_out, &row, s);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_scalar_mul_2_avx2<const N: usize>(row_out: &mut Vector<N>,
                                                   row_in: &Vector<N>,
                                                   s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr = row_in.as_ptr(); // *const u8
        let ptr_out = row_out.as_ptr(); // *const u8
        let b = _mm_set1_epi8(s.0 as i8);

        for col in (0..N).step_by(32) {
            let a1: __m128i = _mm_loadu_si128(ptr.add(col +  0) as *const __m128i);
            let a2: __m128i = _mm_loadu_si128(ptr.add(col + 16) as *const __m128i);
            let t2: __m256i = gf127_mul_avx2(a1, a2, b, b); 

            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t2);
        }
    }
}

/// c[i] = a[i]*b[i] mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row_out = Vector::<N>::new();
/// let row1 = Vector::<N>::from_u8(1);
/// let row2 = Vector::<N>::from_u8(1);
/// gf127_row_mul_avx2(&mut row_out, &row1, &row2);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_mul_avx2<const N: usize>(row_out: &mut Vector<N>,
                                          in1: &Vector<N>,
                                          in2: &Vector<N>) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr1 = in1.as_ptr(); // *const u8
        let ptr2 = in2.as_ptr(); // *const u8
        let ptr_out = row_out.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let a1: __m128i = _mm_loadu_si128(ptr1.add(col +  0) as *const __m128i);
            let a2: __m128i = _mm_loadu_si128(ptr1.add(col + 16) as *const __m128i);
            let b1: __m128i = _mm_loadu_si128(ptr2.add(col +  0) as *const __m128i);
            let b2: __m128i = _mm_loadu_si128(ptr2.add(col + 16) as *const __m128i);
            let t2: __m256i = gf127_mul_avx2(a1, a2, b1, b2); 
            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t2);
        }
    }
}


/// c[i] == 0 for all i  < n
/// # Examples
/// ```
/// const N: usize = 128;
/// let row1 = Vector::<N>::from_u8(1);
/// Vector::<N>::count_zero(&row1);
/// ```
///
/// # Parameters
/// - `a`: row 
///
/// # return 
///     number of zeros in row
#[target_feature(enable = "avx,avx2")]
pub fn gf127_count_zero_avx2<const N: usize>(row: &Vector<N>) -> u32 {
    assert!(N % 32 == 0);
    unsafe {
        let mut zero: __m256i = _mm256_setzero_si256();
        let mut ask: __m256i = _mm256_setzero_si256();
        let mut mask: __m256i = _mm256_set1_epi8(1);
        let ptr1 = row.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let a: __m256i = _mm256_loadu_si256(ptr1.add(col +  0) as *const __m256i);

            let t2: __m256i = gf127_mul_avx2(a1, a2, b1, b2); 
            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t2);
        }

    }
    let mut c: u32 = 0;
    for col in 0..N {
        if row[col] == Fq(0) {
            c += 1;
        }
    }
    return c;
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

    #[test]
    fn test_gf127_mul_avx2() { 
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a1 = _mm_set1_epi8(i as i8);
                    let a2 = _mm_set1_epi8(i as i8);
                    let b1 = _mm_set1_epi8(j as i8);
                    let b2 = _mm_set1_epi8(j as i8);
                    let c = gf127_mul_avx2(a1, a2, b1, b2);
                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, ((i as u16 * j as u16) % 127) as u8);
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_acc_avx2() { 
        const N: usize = 128;
        let mut row = Vector::<N>::from_u8(1);
        unsafe {
            let a1 = gf127_row_acc_avx2(&mut row);
            assert_eq!(a1, (N % 127) as u8);
        }
    }


    #[test]
    fn test_gf127_row_acc_inv_avx2() { 
        const N: usize = 128;
        let mut row = Vector::<N>::from_u8(1);
        unsafe {
            let a1 = gf127_row_acc_inv_avx2(&mut row);
            assert_eq!(a1, (N % 127) as u8);
        }
    }

    #[test]
    fn test_gf127_row_scalar_mul_avx2() { 
        const N: usize = 128;
        let mut row = Vector::<N>::from_u8(1);
        unsafe {
            gf127_row_scalar_mul_avx2(&mut row, Fq(2));
            for i in 0..N {
                assert_eq!(Fq(2), row[i]);
            }
        }
    }
    
    #[test]
    fn test_gf127_row_scalar_mul_2_avx2() { 
        const N: usize = 128;
        let mut row_out = Vector::<N>::from_u8(1);
        let row = Vector::<N>::from_u8(0);
        unsafe {
            gf127_row_scalar_mul_2_avx2(&mut row_out, &row, Fq(2));
            for i in 0..N {
                assert_eq!(Fq(0), row_out[i]);
            }
        }
    }

    #[test]
    fn test_gf127_row_mul_avx2() { 
        const N: usize = 128;
        let mut row_out = Vector::<N>::new();
        let row1 = Vector::<N>::from_u8(1);
        let row2 = Vector::<N>::from_u8(1);

        unsafe {
            gf127_row_mul_avx2(&mut row_out, &row1, &row2);
            for i in 0..N {
                assert_eq!(Fq(1), row_out[i]);
            }
        }
    }
}
