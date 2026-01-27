#![allow(unused_imports)]


use crate::{
    fq::Fq,
    fq::FQ127_INV_TABLE,
    vector::Vector
};

use std::arch::x86_64::{
    __m128i, __m256i, _mm_extract_epi8, _mm_loadu_si128, _mm_set1_epi8, _mm_storeu_si128, _mm256_add_epi8, _mm256_and_si256, _mm256_blendv_epi8, _mm256_cmpeq_epi8, _mm256_cvtepu8_epi16, _mm256_extract_epi8, _mm256_extract_epi32, _mm256_load_si256, _mm256_loadu_si256, _mm256_min_epu8, _mm256_movemask_epi8, _mm256_mullo_epi16, _mm256_or_si256, _mm256_packs_epi16, _mm256_permute2x128_si256, _mm256_permute4x64_epi64, _mm256_set1_epi8, _mm256_set1_epi16, _mm256_setzero_si256, _mm256_srli_epi16, _mm256_srli_epi32, _mm256_srli_epi64, _mm256_srli_si256, _mm256_store_si256, _mm256_storeu_si256, _mm256_sub_epi8, _mm512_setzero
};
use std::arch::x86_64::{
    __mmask64,
    __m512i, 
    _MM_CMPINT_EQ,
    _mm512_set1_epi8, _mm512_set1_epi16,
    _mm512_add_epi8, _mm512_add_epi16,
    _mm512_sub_epi8, _mm512_sub_epi16,
    _mm512_min_epu8,
    _mm512_extracti32x4_epi32,
    _mm512_load_si512, _mm512_store_si512,
    _mm512_permutex2var_epi8,
    _mm512_cvtepu8_epi16, _mm512_cvtepi16_epi8, 
    _mm512_mullo_epi16, _mm512_mulhi_epi16, _mm512_mulhi_epu16,
    _mm512_castsi512_si256, _mm512_castsi256_si512,
    _mm512_extracti32x8_epi32,
    _mm512_setzero_si512,
    _mm512_cmp_epu8_mask, _mm256_cmp_epu8_mask,
};

/// using blend instruction
/// c[i] = a[i] mod q,  a[i] < 2*127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(i as i8);
/// let c = gf127_red_u256(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
#[inline]
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
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx,avx2")]
fn gf127_red_u256_v2(a: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let b = _mm256_sub_epi8(a, q);
    let t3 = _mm256_min_epu8(b, a);
    t3
}

/// c[i] = a[i] mod q,  a[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm512_set1_epi8(i as i8);
/// let c = gf127_red_u16_u512(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
fn gf127_red_u16_u512(a: __m512i) -> __m512i {
    let c01 = _mm512_set1_epi16(0x01);
    let c7f = _mm512_set1_epi16(0x7f);
    let c516 = _mm512_set1_epi16(516);

    let mut tmp: __m512i = _mm512_add_epi16(a, c01);
    tmp = _mm512_mulhi_epu16(tmp, c516);
    tmp = _mm512_mullo_epi16(tmp, c7f);
    tmp = _mm512_sub_epi16(a, tmp);

    return tmp;
}

/// c[i] = a[i] mod q,  a[i] < 127, for i in 0..64
/// # Examples
/// ```
/// let a = _mm512_set1_epi8(i as i8);
/// let c = gf127_red_u512(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
fn gf127_red_u512(a: __m512i) -> __m512i {
    let q = _mm512_set1_epi8(0x7f);
    let t2 = _mm512_sub_epi8(a, q);
    let t3 = _mm512_min_epu8(a, t2);
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
#[inline]
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
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx,avx2")]
fn gf127_add_u256_v2(a: __m256i, b: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_add_epi8(a, b);
    let t2 = _mm256_sub_epi8(t1, q);
    let t3 = _mm256_min_epu8(t1, t2);
    t3
}

/// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm512_set1_epi8(i as i8);
/// let b = _mm512_set1_epi8(j as i8);
/// let c = gf127_add_u512(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
fn gf127_add_u512(a: __m512i, b: __m512i) -> __m512i {
    let q = _mm512_set1_epi8(0x7f);
    let t1 = _mm512_add_epi8(a, b);
    let t2 = _mm512_sub_epi8(t1, q);
    let t3 = _mm512_min_epu8(t1, t2);
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
#[inline]
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
/// let c = gf127_sub_u512a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_sub_u256_v2(a: __m256i, b: __m256i) -> __m256i {
    let q  = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_sub_epi8(a, b);
    let t2 = _mm256_add_epi8(t1, q);
    let t3 = _mm256_min_epu8(t1, t2);
    t3
}

/// c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm512_set1_epi8(i as i8);
/// let b = _mm512_set1_epi8(j as i8);
/// let c = gf127_add_u512(a, b);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
fn gf127_sub_u512(a: __m512i, b: __m512i) -> __m512i {
    let q = _mm512_set1_epi8(0x7f);
    let t1 = _mm512_sub_epi8(a, b);
    let t2 = _mm512_add_epi8(t1, q);
    let t3 = _mm512_min_epu8(t1, t2);
    t3
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
/// - `a1`: first/lower half of vector a 
/// - `a2`: second/higher half of vector a
/// - `b1`: first/lower half of vector b
/// - `b2`: second/higher half of vector b
#[inline]
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

/// c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm256_set1_epi8(0 as i8);
/// let b = _mm256_set1_epi8(0 as i8);
/// let c = gf127_mul_avx512(a, b);
/// ```
///
/// # Parameters
/// - `a`: vector a 
/// - `b`: vector a 
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn gf127_mul_avx512(a: __m256i, b: __m256i) -> __m256i {
    unsafe {
        let a_: __m512i = _mm512_cvtepu8_epi16(a);
        let b_: __m512i = _mm512_cvtepu8_epi16(b);
        return gf127_mul_u16_avx512(a_, b_);
    }
}

/// c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a: __m512i = _mm512_set1_epi16(0 as i8);
/// let b: __m512i = _mm512_set1_epi16(0 as i8);
/// let c: __m256i = gf127_mul_u16_avx512(a, b);
/// ```
///
/// # Parameters
/// - `a`: vector a 
/// - `b`: vector a 
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn gf127_mul_u16_avx512(a: __m512i, b: __m512i) -> __m256i {
    let c7f = _mm512_set1_epi16(127);
    let c01 = _mm512_set1_epi16(1);
    let c516 = _mm512_set1_epi16(516);

    let mut acc = _mm512_mullo_epi16(a, b);
    let mut tmp = _mm512_add_epi16(acc, c01);
    tmp = _mm512_mulhi_epi16(tmp, c516);
    tmp = _mm512_mullo_epi16(tmp, c7f);
    acc = _mm512_sub_epi16(acc, tmp);
    let t = _mm512_cvtepi16_epi8(acc);
    return t;
}

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
#[inline]
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

///  sum_{i=0}^{i < 32} v[i] mod q, v[i] < 127, for i in 0..32
/// # Examples
/// ```
/// let a = _mm512_set1_epi8(1 as i8);
/// let c = gf127_hadd_avx512(a);
/// ```
///
/// # Parameters
/// - `a`: first addend
/// - `b`: second addend
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
fn gf127_hadd_avx512(v: __m512i) -> u8 {
    unsafe {
        let lo = _mm512_castsi512_si256(v);
        let hi = _mm512_extracti32x8_epi32(v, 1);

        let x0 = _mm256_add_epi8(lo, hi);
        let x1 = gf127_red_u256(x0);
        return gf127_hadd_avx2(x1);
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_add_avx2<const N: usize>(c: &mut Vector<N>,
                                          a: &Vector<N>,
                                          b: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr2 = b.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr2.add(col) as *const __m256i);
            let tc = gf127_add_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_add2_avx2<const N: usize>(c: &mut Vector<N>,
                                           a: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr_out.add(col) as *const __m256i);
            let tc = gf127_add_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_add_avx512<const N: usize>(c: &mut Vector<N>,
                                            a: &Vector<N>,
                                            b: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr2 = b.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let ta: __m512i = _mm512_load_si512(ptr1.add(col) as *const __m512i);
            let tb: __m512i = _mm512_load_si512(ptr2.add(col) as *const __m512i);
            let tc = gf127_add_u512(ta, tb); 
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, tc);
            off += 64;
        }

        for col in (off..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr2.add(col) as *const __m256i);
            let tc = gf127_add_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_add2_avx512<const N: usize>(c: &mut Vector<N>,
                                            a: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let ta: __m512i = _mm512_load_si512(ptr1.add(col) as *const __m512i);
            let tb: __m512i = _mm512_load_si512(ptr_out.add(col) as *const __m512i);
            let tc = gf127_add_u512(ta, tb); 
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, tc);
            off += 64;
        }

        for col in (off..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr_out.add(col) as *const __m256i);
            let tc = gf127_add_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}


/// TODO doc
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_sub_avx2<const N: usize>(c: &mut Vector<N>,
                                          a: &Vector<N>,
                                          b: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr2 = b.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr2.add(col) as *const __m256i);
            let tc = gf127_sub_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_sub2_avx2<const N: usize>(c: &mut Vector<N>,
                                           a: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr_out.add(col) as *const __m256i);
            let tc = gf127_sub_u256(tb, ta);
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_sub_avx512<const N: usize>(c: &mut Vector<N>,
                                            a: &Vector<N>,
                                            b: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr2 = b.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let ta: __m512i = _mm512_load_si512(ptr1.add(col) as *const __m512i);
            let tb: __m512i = _mm512_load_si512(ptr2.add(col) as *const __m512i);
            let tc = gf127_sub_u512(ta, tb); 
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, tc);
            off += 64;
        }

        for col in (off..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr2.add(col) as *const __m256i);
            let tc = gf127_sub_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_sub2_avx512<const N: usize>(c: &mut Vector<N>,
                                             a: &Vector<N>) {
    unsafe {
        let ptr1 = a.as_ptr(); // *const u8
        let ptr_out = c.as_ptr(); // *const u8

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let ta: __m512i = _mm512_load_si512(ptr1.add(col) as *const __m512i);
            let tb: __m512i = _mm512_load_si512(ptr_out.add(col) as *const __m512i);
            let tc = gf127_sub_u512(ta, tb); 
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, tc);
            off += 64;
        }

        for col in (off..N).step_by(32) {
            let ta: __m256i = _mm256_load_si256(ptr1.add(col) as *const __m256i);
            let tb: __m256i = _mm256_load_si256(ptr_out.add(col) as *const __m256i);
            let tc = gf127_sub_u256(ta, tb); 
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, tc);
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
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_mul_avx2<const N: usize>(row_out: &mut Vector<N>,
                                          in1: &Vector<N>,
                                          in2: &Vector<N>) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr1 = in1.as_ptr();
        let ptr2 = in2.as_ptr();
        let ptr_out = row_out.as_ptr();

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

/// c[i] = c[i]*a[i] mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row_out = Vector::<N>::new();
/// let row1 = Vector::<N>::from_u8(1);
/// gf127_row_mul2_avx2(&mut row_out, &row1);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_mul2_avx2<const N: usize>(row_out: &mut Vector<N>,
                                           in1: &Vector<N>) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr1 = in1.as_ptr();
        let ptr_out = row_out.as_ptr();

        for col in (0..N).step_by(32) {
            let a1: __m128i = _mm_loadu_si128(ptr1.add(col +  0) as *const __m128i);
            let a2: __m128i = _mm_loadu_si128(ptr1.add(col + 16) as *const __m128i);
            let b1: __m128i = _mm_loadu_si128(ptr_out.add(col +  0) as *const __m128i);
            let b2: __m128i = _mm_loadu_si128(ptr_out.add(col + 16) as *const __m128i);
            let t2: __m256i = gf127_mul_avx2(a1, a2, b1, b2);
            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t2);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_mul_avx512<const N: usize>(row: &mut Vector<N>,
                                            in1: &Vector<N>,
                                            in2: &Vector<N>) {
    unsafe {
        let ptr1 = in1.as_ptr();
        let ptr2 = in2.as_ptr();
        let ptr_out = row.as_ptr();

        for col in (0..N).step_by(32) {
            let a: __m256i = _mm256_loadu_si256(ptr1.add(col) as *const __m256i);
            let b: __m256i = _mm256_loadu_si256(ptr2.add(col) as *const __m256i);
            let t = gf127_mul_avx512(a, b);
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, t);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_mul2_avx512<const N: usize>(row: &mut Vector<N>,
                                             in1: &Vector<N>) {
    unsafe {
        let ptr1 = in1.as_ptr();
        let ptr_out = row.as_ptr();

        for col in (0..N).step_by(32) {
            let a: __m256i = _mm256_loadu_si256(ptr1.add(col) as *const __m256i);
            let b: __m256i = _mm256_loadu_si256(ptr_out.add(col) as *const __m256i);
            let t = gf127_mul_avx512(a, b);
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, t);
        }
    }
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
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_acc_avx2<const N: usize>(row: &Vector<N>) -> u8 {
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

/// \return sum_i..N a[i] mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let a1 = gf127_row_acc_avx512(&mut row);
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_acc_avx512<const N: usize>(row: &Vector<N>) -> u8 {
    assert!(N % 32 == 0);

    unsafe {
        let ptr = row.as_ptr(); // *const u8
        
        let t0: __m256i = _mm256_loadu_si256(ptr as *const __m256i);
        let mut acc = _mm512_cvtepu8_epi16(t0);

        for col in (32..N).step_by(32) {
            let t0: __m256i = _mm256_loadu_si256(ptr.add(col) as *const __m256i);
            let t1 = _mm512_cvtepu8_epi16(t0);
            acc = _mm512_add_epi16(acc, t1);
        }

        return gf127_hadd_avx512(acc);
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
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_acc_inv_avx2<const N: usize>(row: &Vector<N>) -> u8 {
    assert!(N % 32 == 0);
    let mut inv_data: Vector<N> = Vector::new();
    for col in 0..N {
        inv_data[col] = Fq::inv_non_ct(row[col]);
    }

    return gf127_row_acc_avx2(&mut inv_data);
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
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_acc_inv_avx512<const N: usize>(row: &Vector<N>) -> u8 {
    assert!(N % 32 == 0);

    unsafe {
        let ptr_inv = FQ127_INV_TABLE.as_ptr(); // *const u8
        let t1 = _mm512_load_si512(ptr_inv.add( 0) as *const __m512i);
        let t2 = _mm512_load_si512(ptr_inv.add(64) as *const __m512i);

        let mut inv_data: Vector<N> = Vector::new();

        let ptr = row.as_ptr(); // *const u8
        let ptr_out = inv_data.as_ptr(); // *const u8
        
        let mut off = 0;
        for col in (0..N).step_by(64) {
            let a = _mm512_load_si512(ptr.add(col) as *const __m512i);
            let t = _mm512_permutex2var_epi8(t1, a, t2);
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, t);
            off += 64;
        }
        
        for col in (off..N).step_by(64) {
            let b = _mm256_load_si256(ptr.add(col) as *const __m256i);
            let a = _mm512_castsi256_si512(b);
            let t = _mm512_permutex2var_epi8(t1, a, t2);
            let u = _mm512_castsi512_si256(t);

            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, u);
            off += 64;
        }
        
        return gf127_row_acc_avx512(&mut inv_data);
    };
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
#[allow(dead_code)]
#[inline]
#[target_feature(enable = "avx,avx2")]
unsafe fn gf127_row_scalar_mul_avx2_(a: [Fq; 32], b: __m256i) -> __m256i {
    unsafe {
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
}

/// # Examples
/// ```
/// ```
///
/// # Parameters
/// - `s`: input scalar value
///
/// # Returns
/// - [0*s, ..., 127 * s]
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn gf127_scalar_mul_gen_avx512(s: Fq) -> [__m512i; 2]{
    unsafe {
        let ptr = FQ127_INV_TABLE.as_ptr(); // *const u8
        let ret: [__m512i; 2] = [
            _mm512_load_si512(ptr.add( 0 + 128 * s.0 as usize) as *const __m512i),
            _mm512_load_si512(ptr.add(64 + 128 * s.0 as usize) as *const __m512i)
        ];

        return ret;
    }
}

/// # Examples
/// ```
/// ```
///
/// # Parameters
/// - `s`:
///
/// # Returns
/// - [0*s, ..., 127 * s]
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn gf127_scalar_mul_avx512(a: __m512i,
                                  table: [__m512i; 2]) -> __m512i {
    unsafe {
        return _mm512_permutex2var_epi8(table[0], a, table[1]);
    }
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
#[inline]
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

/// \return a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_avx512(&mut row), s;
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_scalar_mul_avx512<const N: usize>(row: &mut Vector<N>,
                                                s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr = row.as_ptr(); // *const u8
        let b = _mm512_set1_epi8(s.0 as i8);

        for col in (0..N).step_by(32) {
            let a0: __m256i = _mm256_load_si256(ptr.add(col +  0) as *const __m256i);
            let a1 = _mm512_cvtepu8_epi16(a0);
            let t0: __m256i = gf127_mul_u16_avx512(a1, b); 
            _mm256_storeu_si256(ptr.add(col) as *mut __m256i, t0);
        }
    }
}

/// \return a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_avx512(&mut row), s;
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_scalar_mul_avx512_non_ct<const N: usize>(row: &mut Vector<N>,
                                                          s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let tab = gf127_scalar_mul_gen_avx512(s);
        let ptr = row.as_ptr(); // *const u8

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let a0: __m512i = _mm512_load_si512(ptr.add(col +  0) as *const __m512i);
            let a1 = gf127_scalar_mul_avx512(a0, tab);
            _mm512_store_si512(ptr.add(col) as *mut __m512i, a1);
            off += 64
        }
        
        for col in (off..N).step_by(32) {
            let a0: __m256i = _mm256_load_si256(ptr.add(col +  0) as *const __m256i);
            let a1 = _mm512_castsi256_si512(a0);
            let a2 = gf127_scalar_mul_avx512(a1, tab);
            let a3 = _mm512_castsi512_si256(a2);
            _mm256_store_si256(ptr.add(col) as *mut __m256i, a3);
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
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_scalar_mul_2_avx2<const N: usize>(row_out: &mut Vector<N>,
                                                   row_in: &Vector<N>,
                                                   s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr = row_in.as_ptr();
        let ptr_out = row_out.as_ptr();
        let b = _mm_set1_epi8(s.0 as i8);

        for col in (0..N).step_by(32) {
            let a1: __m128i = _mm_loadu_si128(ptr.add(col +  0) as *const __m128i);
            let a2: __m128i = _mm_loadu_si128(ptr.add(col + 16) as *const __m128i);
            let t2: __m256i = gf127_mul_avx2(a1, a2, b, b); 

            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t2);
        }
    }
}

/// TODO test
/// \return c[i] = a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_avx512(&mut row), s;
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_scalar_mul_2_avx512<const N: usize>(row_out: &mut Vector<N>,
                                                     row_in: &Vector<N>,
                                                     s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let ptr_out = row_out.as_ptr(); 
        let ptr_in = row_in.as_ptr(); // *const u8
        let b = _mm512_set1_epi8(s.0 as i8);

        for col in (0..N).step_by(32) {
            let a0: __m256i = _mm256_load_si256(ptr_in.add(col +  0) as *const __m256i);
            let a1 = _mm512_cvtepu8_epi16(a0);
            let t0: __m256i = gf127_mul_u16_avx512(a1, b); 
            _mm256_storeu_si256(ptr_out.add(col) as *mut __m256i, t0);
        }
    }
}

/// TODO test
/// \return a[i]*s mod q, a[i] < 127 for i in 0..N
/// # Examples
/// ```
/// const N: usize = 128;
/// let mut row = Vector::<N>::from_u8(1);
/// let s: u7 = 2;
/// gf127_row_scalar_mul_avx512(&mut row), s;
/// ```
///
/// # Parameters
/// - `a`: row to accumulate
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_scalar_mul_2_avx512_non_ct<const N: usize>(row_out: &mut Vector<N>,
                                                            row_in: &Vector<N>,
                                                            s: Fq) {
    assert!(N % 32 == 0);
    unsafe {
        let tab = gf127_scalar_mul_gen_avx512(s);
        let ptr_out = row_out.as_ptr();
        let ptr_in = row_in.as_ptr();

        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let a0: __m512i = _mm512_load_si512(ptr_in.add(col +  0) as *const __m512i);
            let a1 = gf127_scalar_mul_avx512(a0, tab);
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, a1);
            off += 64
        }
        
        for col in (off..N).step_by(32) {
            let a0: __m256i = _mm256_load_si256(ptr_in.add(col +  0) as *const __m256i);
            let a1 = _mm512_castsi256_si512(a0);
            let a2 = gf127_scalar_mul_avx512(a1, tab);
            let a3 = _mm512_castsi512_si256(a2);
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, a3);
        }
    }
}

/// TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_inv_avx512<const N: usize>(row_out: &mut Vector<N>,
                                            row: &Vector<N>) {
    unsafe {
        let tab = FQ127_INV_TABLE.as_ptr();
        let ptr = row.as_ptr();
        let ptr_out = row_out.as_ptr();

        let t0 = _mm512_load_si512(tab.add( 0) as *const __m512i);
        let t1 = _mm512_load_si512(tab.add(64) as *const __m512i);
       
        let mut off: usize = 0;
        for col in (0..N).step_by(64) {
            let a: __m512i = _mm512_load_si512(ptr.add(col) as *const __m512i);
            let t = _mm512_permutex2var_epi8(t0, a, t1);
            _mm512_store_si512(ptr_out.add(col) as *mut __m512i, t);
            off += 64;
        }

        for col in (off..N).step_by(64) {
            let a: __m256i = _mm256_load_si256(ptr.add(col) as *const __m256i);
            let b: __m512i = _mm512_castsi256_si512(a);
            let c = _mm512_permutex2var_epi8(t0, b, t1);
            let t: __m256i = _mm512_castsi512_si256(c);
            _mm256_store_si256(ptr_out.add(col) as *mut __m256i, t);
            off += 64;
        }
    }
}

/// c[i] == 0 for all i  < n
/// # Examples
/// ```
/// const N: usize = 128;
/// let row1 = Vector::<N>::from_u8(1);
/// gf127_row_contains_zero_avx2(&row1);
/// ```
///
/// # Parameters
/// - `a`: row 
///
/// # return 
///     1 if row contains a zer 
///     0 else
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_contains_zero_avx2<const N: usize>(row: &Vector<N>) -> bool {
    assert!(N % 32 == 0);
    unsafe {
        let zero: __m256i = _mm256_setzero_si256();
        let mut acc: __m256i = _mm256_setzero_si256();
        let ptr1 = row.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let a0: __m256i = _mm256_loadu_si256(ptr1.add(col +  0) as *const __m256i);
            let a1: __m256i = _mm256_cmpeq_epi8(a0, zero);
            acc = _mm256_or_si256(acc, a1);
        }

        let t: i32 = _mm256_movemask_epi8(acc);
        return t != 0;
    }
}

///  TODO doc
#[inline]
#[target_feature(enable = "avx512f,avx512bw")]
pub fn gf127_row_contains_zero_avx512<const N: usize>(row: &Vector<N>) -> bool {
    assert!(N % 32 == 0);
    unsafe {
        let ptr1 = row.as_ptr(); // *const u8
        let zero = _mm512_setzero_si512();

        let mut off: usize = 0;
        let mut acc: __mmask64 = 0;
        for col in (0..N).step_by(64) {
            let t0: __m512i = _mm512_load_si512(ptr1.add(col) as *const __m512i);
            acc |= _mm512_cmp_epu8_mask(t0, zero, _MM_CMPINT_EQ);
            off += 64;
        }
        for col in (off..N).step_by(32) {
            let t0: __m256i = _mm256_loadu_si256(ptr1.add(col) as *const __m256i);
            acc |= _mm256_cmp_epu8_mask(t0, _mm512_castsi512_si256(zero), _MM_CMPINT_EQ) as u64;
        }
        return acc != 0;
    }
}

/// c[i] == 0 for all i  < n
/// # Examples
/// ```
/// const N: usize = 128;
/// let row1 = Vector::<N>::from_u8(1);
/// gf127_row_count_zero_avx2(&row1);
/// ```
///
/// # Parameters
/// - `a`: row 
///
/// # return 
///     number of zeros in row
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_row_count_zero_avx2<const N: usize>(row: &Vector<N>) -> u32 {
    assert!(N % 32 == 0);
    unsafe {
        let zero: __m256i = _mm256_setzero_si256();
        let mask: __m256i = _mm256_set1_epi8(1);
        let mut acc: __m256i = _mm256_setzero_si256();
        let ptr1 = row.as_ptr(); // *const u8

        for col in (0..N).step_by(32) {
            let a0: __m256i = _mm256_loadu_si256(ptr1.add(col +  0) as *const __m256i);
            let a1: __m256i = _mm256_cmpeq_epi8(a0, zero);
            let a2: __m256i = _mm256_and_si256(a1, mask);
            acc = _mm256_add_epi8(acc, a2);
        }

        return gf127_hadd_avx2(acc) as u32;
    }
}


/// TODO
#[inline]
#[target_feature(enable = "avx,avx2")]
pub fn gf127_histogram_avx2<const N: usize>(out: &[u8; 128], row: &Vector<N>) -> u32 {
    return 0;
}




#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 128;

    #[test]
    fn test_gf127_red_u256() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        unsafe {
            for i in 0..254u8 {
                let a = _mm256_set1_epi8(i as i8);
                let c = gf127_red_u256(a);

                let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                assert_eq!(t, (i % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_red_u256_v2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        unsafe {
            for i in 0..254u8 {
                let a = _mm256_set1_epi8(i as i8);
                let c = gf127_red_u256_v2(a);

                let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                assert_eq!(t, (i % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_red_u512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..16129i16 {
                let a = _mm512_set1_epi16(i);
                let c = gf127_red_u16_u512(a);
                let t: u8 = (_mm_extract_epi8::<0>(_mm512_extracti32x4_epi32::<0>(c))) as u8;
                assert_eq!(t, (i % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_red_u16_u512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..254u8 {
                let a = _mm512_set1_epi8(i as i8);
                let c = gf127_red_u512(a);

                let t: u8 = (_mm_extract_epi8::<0>(_mm512_extracti32x4_epi32::<0>(c))) as u8;
                assert_eq!(t, (i % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_add_u256() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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

    // #[target_feature(enable = "avx512f,avx512bw")]
    #[test]
    fn test_gf127_add_u512() { 
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm512_set1_epi8(i as i8);
                    let b = _mm512_set1_epi8(j as i8);
                    let c = gf127_add_u512(a, b);

                    let t: u8 = (_mm_extract_epi8::<0>(_mm512_extracti32x4_epi32::<0>(c))) as u8;
                    assert_eq!(t, (i+j) % 127);
                }
            }
        }
    }

    #[test]
    fn test_gf127_sub_u256() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
    fn test_gf127_sub_u512() { 
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm512_set1_epi8(i as i8);
                    let b = _mm512_set1_epi8(j as i8);
                    let c = gf127_sub_u512(a, b);

                    let t: u8 = (_mm_extract_epi8::<0>(_mm512_extracti32x4_epi32::<0>(c))) as u8;
                    assert_eq!(t, (i+127-j) % 127);
                }
            }
        }
    }

    #[test]
    fn test_gf127_mul_avx2() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
    fn test_gf127_mul_avx512() { 
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..127u8 {
                for j in 0..127u8 {
                    let a = _mm256_set1_epi8(i as i8);
                    let b = _mm256_set1_epi8(j as i8);
                    let c = gf127_mul_avx512(a, b);
                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, ((i as u16 * j as u16) % 127) as u8);
                }
            }
        }
    }

    #[test]
    fn test_gf127_mul_u16_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..127i16 {
                for j in 0..127i16 {
                    let a = _mm256_set1_epi16(i);
                    let b = _mm256_set1_epi16(j);
                    let c = gf127_mul_avx512(a, b);
                    let t: u8 = (_mm256_extract_epi32::<0>(c) & 0xFF) as u8;
                    assert_eq!(t, ((i * j) % 127) as u8);
                }
            }
        }
    }

    #[test]
    fn test_gf127_hadd_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        unsafe {
            for i in 0..127u8 {
                let a1 = _mm256_set1_epi8(i as i8);
                let c = gf127_hadd_avx2(a1);
                assert_eq!(c, ((i as u16 * 32u16) % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_hadd_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        unsafe {
            for i in 0..127u8 {
                let a1 = _mm512_set1_epi8(i as i8);
                let c = gf127_hadd_avx512(a1);
                assert_eq!(c, ((i as u16 * 64u16) % 127) as u8);
            }
        }
    }

    #[test]
    fn test_gf127_row_add_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut c  = Vector::<N>::from_u8(1);
        for i in 0..127u8 {
            for j in 0..127u8 {
                let a = Vector::<N>::from_u8(i);
                let b = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_add_avx2(&mut c, &a, &b);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_add2_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        for i in 0..127u8 {
            for j in 0..127u8 {
                let mut c  = Vector::<N>::from_u8(i);
                let a = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_add2_avx2(&mut c, &a);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_add_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        let mut c  = Vector::<N>::from_u8(1);
        for i in 0..127u8 {
            for j in 0..127u8 {
                let a = Vector::<N>::from_u8(i);
                let b = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_add_avx512(&mut c, &a, &b);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_add2_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        for i in 0..127u8 {
            for j in 0..127u8 {
                let mut c  = Vector::<N>::from_u8(i);
                let a = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_add2_avx512(&mut c, &a);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+j) % 127);
                    }
                }
            }
        }
    }


    #[test]
    fn test_gf127_row_sub_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut c  = Vector::<N>::from_u8(1);
        for i in 0..127u8 {
            for j in 0..127u8 {
                let a = Vector::<N>::from_u8(i);
                let b = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_sub_avx2(&mut c, &a, &b);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+127-j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_sub2_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        for i in 0..127u8 {
            for j in 0..127u8 {
                let mut c  = Vector::<N>::from_u8(i);
                let a = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_sub2_avx2(&mut c, &a);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+127-j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_sub_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        let mut c  = Vector::<N>::from_u8(1);
        for i in 0..127u8 {
            for j in 0..127u8 {
                let a = Vector::<N>::from_u8(i);
                let b = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_sub_avx512(&mut c, &a, &b);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+127-j) % 127);
                    }
                }
            }
        }
    }

    #[test]
    fn test_gf127_row_sub2_avx512() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        for i in 0..127u8 {
            for j in 0..127u8 {
                let mut c  = Vector::<N>::from_u8(i);
                let a = Vector::<N>::from_u8(j);
                unsafe {
                    gf127_row_sub2_avx512(&mut c, &a);
                    for k in 0..N {
                        assert_eq!(c[k].0, (i+127-j) % 127);
                    }
                }
            }
        }
    }


    #[test]
    fn test_gf127_row_mul_avx2() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }

        for i in 0..127u8 {
            for j in 0..127u8 {
                let mut row_out = Vector::<N>::new();
                let row1 = Vector::<N>::from_u8(i);
                let row2 = Vector::<N>::from_u8(j);

                unsafe {
                    gf127_row_mul_avx2(&mut row_out, &row1, &row2);
                    for i in 0..N {
                        assert_eq!(row_out[i].0, (((i as u16) * (j as u16)) % 127) as u8);
                    }
                }
            }
        }
    }






    #[test]
    fn test_gf127_row_scalar_mul_avx2() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
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
    fn test_gf127_row_contains_zero_avx2() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut row_out = Vector::<N>::from_u8(1);

        unsafe {
            let t1 = gf127_row_contains_zero_avx2(&mut row_out,);
            assert_eq!(false, t1);
        }
       
        row_out[17] = Fq(0);
        unsafe {
            let t1 = gf127_row_contains_zero_avx2(&mut row_out,);
            assert_eq!(true, t1);
        }
    }

    #[test]
    fn test_gf127_row_count_zero_avx2() { 
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut row_out = Vector::<N>::from_u8(1);

        unsafe {
            let t1 = gf127_row_count_zero_avx2(&mut row_out,);
            assert_eq!(0, t1);
        }
       
        row_out[17] = Fq(0);
        unsafe {
            let t1 = gf127_row_count_zero_avx2(&mut row_out,);
            assert_eq!(1, t1);
        }
    }
}
