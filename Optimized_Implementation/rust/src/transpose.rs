//use core::arch::x86_64::*;
//use core::ptr;
//
///// Transpose an 8×8 byte matrix.
/////
///// # Safety
/////
///// - `src` and `dst` must each contain at least 8 rows.
///// - Each row must contain at least 8 bytes.
///// - `src_stride`/`dst_stride` are measured in bytes.
//pub unsafe fn matrix_transpose8x8(
//    dst: *mut u8,
//    src: *const u8,
//    src_stride: usize,
//    dst_stride: usize,
//) {
//    // Load rows.
//    let a0 = ptr::read_unaligned(src.add(0 * src_stride) as *const u64);
//    let a1 = ptr::read_unaligned(src.add(1 * src_stride) as *const u64);
//    let a2 = ptr::read_unaligned(src.add(2 * src_stride) as *const u64);
//    let a3 = ptr::read_unaligned(src.add(3 * src_stride) as *const u64);
//    let a4 = ptr::read_unaligned(src.add(4 * src_stride) as *const u64);
//    let a5 = ptr::read_unaligned(src.add(5 * src_stride) as *const u64);
//    let a6 = ptr::read_unaligned(src.add(6 * src_stride) as *const u64);
//    let a7 = ptr::read_unaligned(src.add(7 * src_stride) as *const u64);
//
//    // 2×2 blocks.
//    let b0 = (a0 & 0x00ff00ff00ff00ff) | ((a1 << 8) & 0xff00ff00ff00ff00);
//    let b1 = (a1 & 0xff00ff00ff00ff00) | ((a0 >> 8) & 0x00ff00ff00ff00ff);
//    let b2 = (a2 & 0x00ff00ff00ff00ff) | ((a3 << 8) & 0xff00ff00ff00ff00);
//    let b3 = (a3 & 0xff00ff00ff00ff00) | ((a2 >> 8) & 0x00ff00ff00ff00ff);
//    let b4 = (a4 & 0x00ff00ff00ff00ff) | ((a5 << 8) & 0xff00ff00ff00ff00);
//    let b5 = (a5 & 0xff00ff00ff00ff00) | ((a4 >> 8) & 0x00ff00ff00ff00ff);
//    let b6 = (a6 & 0x00ff00ff00ff00ff) | ((a7 << 8) & 0xff00ff00ff00ff00);
//    let b7 = (a7 & 0xff00ff00ff00ff00) | ((a6 >> 8) & 0x00ff00ff00ff00ff);
//
//    // 4×4 blocks.
//    let c0 = (b0 & 0x0000ffff0000ffff) | ((b2 << 16) & 0xffff0000ffff0000);
//    let c1 = (b1 & 0x0000ffff0000ffff) | ((b3 << 16) & 0xffff0000ffff0000);
//    let c2 = (b2 & 0xffff0000ffff0000) | ((b0 >> 16) & 0x0000ffff0000ffff);
//    let c3 = (b3 & 0xffff0000ffff0000) | ((b1 >> 16) & 0x0000ffff0000ffff);
//    let c4 = (b4 & 0x0000ffff0000ffff) | ((b6 << 16) & 0xffff0000ffff0000);
//    let c5 = (b5 & 0x0000ffff0000ffff) | ((b7 << 16) & 0xffff0000ffff0000);
//    let c6 = (b6 & 0xffff0000ffff0000) | ((b4 >> 16) & 0x0000ffff0000ffff);
//    let c7 = (b7 & 0xffff0000ffff0000) | ((b5 >> 16) & 0x0000ffff0000ffff);
//
//    // 8×8 block.
//    let d0 = (c0 & 0x00000000ffffffff) | ((c4 << 32) & 0xffffffff00000000);
//    let d1 = (c1 & 0x00000000ffffffff) | ((c5 << 32) & 0xffffffff00000000);
//    let d2 = (c2 & 0x00000000ffffffff) | ((c6 << 32) & 0xffffffff00000000);
//    let d3 = (c3 & 0x00000000ffffffff) | ((c7 << 32) & 0xffffffff00000000);
//    let d4 = (c4 & 0xffffffff00000000) | ((c0 >> 32) & 0x00000000ffffffff);
//    let d5 = (c5 & 0xffffffff00000000) | ((c1 >> 32) & 0x00000000ffffffff);
//    let d6 = (c6 & 0xffffffff00000000) | ((c2 >> 32) & 0x00000000ffffffff);
//    let d7 = (c7 & 0xffffffff00000000) | ((c3 >> 32) & 0x00000000ffffffff);
//
//    // Store rows.
//    ptr::write_unaligned(dst.add(0 * dst_stride) as *mut u64, d0);
//    ptr::write_unaligned(dst.add(1 * dst_stride) as *mut u64, d1);
//    ptr::write_unaligned(dst.add(2 * dst_stride) as *mut u64, d2);
//    ptr::write_unaligned(dst.add(3 * dst_stride) as *mut u64, d3);
//    ptr::write_unaligned(dst.add(4 * dst_stride) as *mut u64, d4);
//    ptr::write_unaligned(dst.add(5 * dst_stride) as *mut u64, d5);
//    ptr::write_unaligned(dst.add(6 * dst_stride) as *mut u64, d6);
//    ptr::write_unaligned(dst.add(7 * dst_stride) as *mut u64, d7);
//}
//
//
//
//const MATRIX_TRANSPOSE_TABLE: [u32; 16] = [
//    0, 8, 4, 12,
//    2, 10, 6, 14,
//    1, 9, 5, 13,
//    3, 11, 7, 15,
//];
//
///// Transpose a 32×32 byte matrix.
/////
///// # Safety
/////
///// - `src_origin` must point to at least 32 rows separated by `src_stride`.
///// - `dst_origin` must point to at least 32 rows separated by `dst_stride`.
///// - CPU must support AVX2.
//#[target_feature(enable = "avx2")]
//pub unsafe fn matrix_transpose_32x32(
//    dst_origin: *mut u8,
//    src_origin: *const u8,
//    src_stride: usize,
//    dst_stride: usize,
//) {
//    let mut t: [__m256i; 32] = core::mem::zeroed();
//
//    // Load
//    for i in 0..32 {
//        t[i] = _mm256_loadu_si256(
//            src_origin.add(i * src_stride) as *const __m256i
//        );
//    }
//
//    // unpack bytes
//    for i in (0..32).step_by(2) {
//        let t0 = _mm256_unpacklo_epi8(t[i], t[i + 1]);
//        let t1 = _mm256_unpackhi_epi8(t[i], t[i + 1]);
//        t[i] = t0;
//        t[i + 1] = t1;
//    }
//
//    // unpack u16
//    for i in (0..32).step_by(4) {
//        let t0 = _mm256_unpacklo_epi16(t[i], t[i + 2]);
//        let t1 = _mm256_unpacklo_epi16(t[i + 1], t[i + 3]);
//        let t2 = _mm256_unpackhi_epi16(t[i], t[i + 2]);
//        let t3 = _mm256_unpackhi_epi16(t[i + 1], t[i + 3]);
//
//        t[i] = t0;
//        t[i + 1] = t1;
//        t[i + 2] = t2;
//        t[i + 3] = t3;
//    }
//
//    // unpack u32
//    for i in (0..32).step_by(8) {
//        let t0 = _mm256_unpacklo_epi32(t[i], t[i + 4]);
//        let t1 = _mm256_unpacklo_epi32(t[i + 1], t[i + 5]);
//        let t2 = _mm256_unpacklo_epi32(t[i + 2], t[i + 6]);
//        let t3 = _mm256_unpacklo_epi32(t[i + 3], t[i + 7]);
//
//        let t4 = _mm256_unpackhi_epi32(t[i], t[i + 4]);
//        let t5 = _mm256_unpackhi_epi32(t[i + 1], t[i + 5]);
//        let t6 = _mm256_unpackhi_epi32(t[i + 2], t[i + 6]);
//        let t7 = _mm256_unpackhi_epi32(t[i + 3], t[i + 7]);
//
//        t[i] = t0;
//        t[i + 1] = t1;
//        t[i + 2] = t2;
//        t[i + 3] = t3;
//        t[i + 4] = t4;
//        t[i + 5] = t5;
//        t[i + 6] = t6;
//        t[i + 7] = t7;
//    }
//
//    // unpack u64
//    for i in 0..8 {
//        let t0 = _mm256_unpacklo_epi64(t[i], t[i + 8]);
//        let t1 = _mm256_unpackhi_epi64(t[i], t[i + 8]);
//        let t2 = _mm256_unpacklo_epi64(t[i + 16], t[i + 24]);
//        let t3 = _mm256_unpackhi_epi64(t[i + 16], t[i + 24]);
//
//        t[i] = t0;
//        t[i + 8] = t1;
//        t[i + 16] = t2;
//        t[i + 24] = t3;
//    }
//
//    // swap 128-bit lanes
//    for i in 0..16 {
//        let t0 = _mm256_permute2x128_si256(t[i], t[i + 16], 0x20);
//        let t1 = _mm256_permute2x128_si256(t[i], t[i + 16], 0x31);
//
//        t[i] = t0;
//        t[i + 16] = t1;
//    }
//
//    // Store first half
//    for i in 0..16 {
//        let pos = MATRIX_TRANSPOSE_TABLE[i] as usize;
//        _mm256_storeu_si256(
//            dst_origin.add(i * dst_stride) as *mut __m256i,
//            t[pos],
//        );
//    }
//
//    // Store second half
//    for i in 0..16 {
//        let pos = MATRIX_TRANSPOSE_TABLE[i] as usize;
//        _mm256_storeu_si256(
//            dst_origin.add((i + 16) * dst_stride) as *mut __m256i,
//            t[pos + 16],
//        );
//    }
//}
//
//
///// Transpose a 64×64 byte matrix.
/////
///// # Safety
/////
///// - CPU must support AVX-512F + AVX-512BW.
///// - `src_origin` and `dst_origin` must reference sufficiently large buffers.
//#[target_feature(enable = "avx512f,avx512bw")]
//pub unsafe fn matrix_transpose_64x64(
//    dst_origin: *mut u8,
//    src_origin: *const u8,
//    src_stride: usize,
//    dst_stride: usize,
//) {
//    let m1 = _mm512_setr_epi64(0, 1, 8, 9, 4, 5, 12, 13);
//    let m2 = _mm512_setr_epi64(2, 3, 10, 11, 6, 7, 14, 15);
//    let m3 = _mm512_setr_epi64(0, 1, 2, 3, 8, 9, 10, 11);
//    let m4 = _mm512_setr_epi64(4, 5, 6, 7, 12, 13, 14, 15);
//
//    let mut t: [__m512i; 64] = core::mem::zeroed();
//
//    // Load
//    for i in 0..64 {
//        t[i] = _mm512_loadu_si512(
//            src_origin.add(i * src_stride) as *const _
//        );
//    }
//
//    // unpack bytes
//    for i in (0..64).step_by(2) {
//        let a = _mm512_unpacklo_epi8(t[i], t[i + 1]);
//        let b = _mm512_unpackhi_epi8(t[i], t[i + 1]);
//        t[i] = a;
//        t[i + 1] = b;
//    }
//
//    // unpack u16
//    for i in (0..64).step_by(4) {
//        let t0 = _mm512_unpacklo_epi16(t[i], t[i + 2]);
//        let t1 = _mm512_unpacklo_epi16(t[i + 1], t[i + 3]);
//        let t2 = _mm512_unpackhi_epi16(t[i], t[i + 2]);
//        let t3 = _mm512_unpackhi_epi16(t[i + 1], t[i + 3]);
//
//        t[i] = t0;
//        t[i + 1] = t1;
//        t[i + 2] = t2;
//        t[i + 3] = t3;
//    }
//
//    // unpack u32
//    for i in (0..64).step_by(8) {
//        let t0 = _mm512_unpacklo_epi32(t[i], t[i + 4]);
//        let t1 = _mm512_unpacklo_epi32(t[i + 1], t[i + 5]);
//        let t2 = _mm512_unpacklo_epi32(t[i + 2], t[i + 6]);
//        let t3 = _mm512_unpacklo_epi32(t[i + 3], t[i + 7]);
//
//        let t4 = _mm512_unpackhi_epi32(t[i], t[i + 4]);
//        let t5 = _mm512_unpackhi_epi32(t[i + 1], t[i + 5]);
//        let t6 = _mm512_unpackhi_epi32(t[i + 2], t[i + 6]);
//        let t7 = _mm512_unpackhi_epi32(t[i + 3], t[i + 7]);
//
//        t[i] = t0;
//        t[i + 1] = t1;
//        t[i + 2] = t2;
//        t[i + 3] = t3;
//        t[i + 4] = t4;
//        t[i + 5] = t5;
//        t[i + 6] = t6;
//        t[i + 7] = t7;
//    }
//
//    // unpack u64
//    for i in 0..8 {
//        let t0 = _mm512_unpacklo_epi64(t[i], t[i + 8]);
//        let t1 = _mm512_unpackhi_epi64(t[i], t[i + 8]);
//        let t2 = _mm512_unpacklo_epi64(t[i + 16], t[i + 24]);
//        let t3 = _mm512_unpackhi_epi64(t[i + 16], t[i + 24]);
//        let t4 = _mm512_unpacklo_epi64(t[i + 32], t[i + 40]);
//        let t5 = _mm512_unpackhi_epi64(t[i + 32], t[i + 40]);
//        let t6 = _mm512_unpacklo_epi64(t[i + 48], t[i + 56]);
//        let t7 = _mm512_unpackhi_epi64(t[i + 48], t[i + 56]);
//
//        t[i] = t0;
//        t[i + 8] = t1;
//        t[i + 16] = t2;
//        t[i + 24] = t3;
//        t[i + 32] = t4;
//        t[i + 40] = t5;
//        t[i + 48] = t6;
//        t[i + 56] = t7;
//    }
//
//    // swap 128-bit limbs
//    for i in 0..16 {
//        let t0 = _mm512_permutex2var_epi64(t[i], m1, t[i + 16]);
//        let t1 = _mm512_permutex2var_epi64(t[i], m2, t[i + 16]);
//        let t2 = _mm512_permutex2var_epi64(t[i + 32], m1, t[i + 48]);
//        let t3 = _mm512_permutex2var_epi64(t[i + 32], m2, t[i + 48]);
//
//        t[i] = t0;
//        t[i + 16] = t1;
//        t[i + 32] = t2;
//        t[i + 48] = t3;
//    }
//
//    // swap 256-bit limbs
//    for i in 0..32 {
//        let t0 = _mm512_permutex2var_epi64(t[i], m3, t[i + 32]);
//        let t1 = _mm512_permutex2var_epi64(t[i], m4, t[i + 32]);
//
//        t[i] = t0;
//        t[i + 32] = t1;
//    }
//
//    // Store
//    for j in 0..4 {
//        let off = j * 16;
//
//        for i in 0..16 {
//            let pos = MATRIX_TRANSPOSE_TABLE[i] as usize;
//
//            _mm512_storeu_si512(
//                dst_origin.add((off + i) * dst_stride) as *mut _,
//                t[off + pos],
//            );
//        }
//    }
//}