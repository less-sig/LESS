use std::arch::x86_64::{
    __m256i,
    _mm256_set1_epi8, 
    _mm256_add_epi8, _mm256_sub_epi8, _mm256_min_epi8
};

#[target_feature(enable = "avx,avx2")]
fn gf127_add_u256(a: __m256i, b: __m256i) -> __m256i {
    let q = _mm256_set1_epi8(0x7f);
    let t1 = _mm256_add_epi8(a, b);
    let t2 = _mm256_sub_epi8(t1, q);
    let t3 = _mm256_min_epi8(t1, t2);
    t3
}

// #[cfg(target_feature = "avx2")]
// unsafe fn gf127_sub_u256(a: __m256i, b: __m256i) -> __m256i {
//     let q  = _mm256_set1_epi8(0x7f);
//     let t1 = _mm256_sub_epi8(a, b);
//     let t2 = _mm256_add_epi8(t1, q);
//     let t3 = _mm256_min_epi8(t1, t2);
//     t3
// } 
// 
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
            let a = _mm256_set1_epi8(0);
            let b = _mm256_set1_epi8(1);
            let c = gf127_add_u256(a, b);
            for i in 0..32 {
                c[i] = 1;
            }
        }
    }
}
