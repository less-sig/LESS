#include <unity.h>

#include "bigint.h"
#include "stm32f4xx_hal.h"
#include "utils.h"

void setUp(void) {
    // set stuff up here
}

void tearDown(void) {
    // clean stuff up here
}

void callback_1ms(void) {
    if (HAL_GetTick() % 512 == 0) {
        LED_toggle();
    }
}

void test_bigint_add_10(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x060C3F3F, 0x025B5966, 0x011064D1, 0x02F0061A, 0x0429069E,
        0x04B4F194, 0x02E1A2CB, 0x01977619, 0x050BB491, 0x002AB375,
    };

    bigint_word_t vc_actual[10];

    bigint_add(vc_actual, va, vb, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_add_20(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100, 0x01A1DFA7, 0x00B06387,
        0x02237B40, 0x0234E302, 0x01D701A2, 0x023BB390, 0x031356EF, 0x03000194, 0x0385BD8F,
        0x036C0259, 0x0131311E, 0x026313F4, 0x02042DAD, 0x009F8A18, 0x0000F995,
    };
    const bigint_word_t vb[] = {
        0x02C1376C, 0x03861208, 0x03F33DDC, 0x016DEE26, 0x03B27CC9, 0x015D5EE3, 0x005CAD9F,
        0x00A8B222, 0x00C710ED, 0x0122BA7C, 0x01D1533B, 0x0053816D, 0x03A74A93, 0x0227BCFE,
        0x0185D94E, 0x03C4B07E, 0x021CEE7F, 0x02CB3647, 0x01ED2C20, 0x00000D51,
    };
    const bigint_word_t vc_expected[] = {
        0x05123DA4, 0x04ABFC76, 0x05038961, 0x04021B44, 0x051B5DC9, 0x02FF3E8A, 0x010D1126,
        0x02CC2D62, 0x02FBF3EF, 0x02F9BC1E, 0x040D06CB, 0x0366D85C, 0x06A74C27, 0x05AD7A8D,
        0x04F1DBA7, 0x04F5E19C, 0x04800273, 0x04CF63F4, 0x028CB638, 0x000106E6,
    };

    bigint_word_t vc_actual[20];

    bigint_add(vc_actual, va, vb, 20);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_sub_positive_result_10(void) {
    const bigint_word_t va[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vc_expected[] = {
        0x016A32CF, 0x000F848A, 0xFEEFCDC7, 0xFDC7ABDE, 0x0157449E,
        0x01713246, 0x0180DBBD, 0xFF73F799, 0x02D6D18F, 0x0013B1D3,
    };

    bigint_word_t vc_actual[10];

    bigint_sub(vc_actual, va, vb, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_sub_negative_result_10(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0xFE95CD31, 0xFFF07B76, 0x01103239, 0x02385422, 0xFEA8BB62,
        0xFE8ECDBA, 0xFE7F2443, 0x008C0867, 0xFD292E71, 0xFFEC4E2D,
    };

    bigint_word_t vc_actual[10];

    bigint_sub(vc_actual, va, vb, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_sub_positive_result_20(void) {
    const bigint_word_t va[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vb[] = {
        0x01C17BF0, 0x01EC012D, 0x03992F9E, 0x01F11716, 0x02F0F7A3, 0x00FCEBFB, 0x03CEAE1A,
        0x02DB42EB, 0x021AE208, 0x01BA7181, 0x0042367E, 0x028BABFA, 0x0381C241, 0x0021E884,
        0x01D70AA5, 0x0352A584, 0x03423231, 0x00126B61, 0x006099A2, 0x0000073A,
    };
    const bigint_word_t vc_expected[] = {
        0xFFC02798, 0xFE952DA8, 0xFF7B8D5A, 0x003DB872, 0xFD720A1E, 0x001799E5, 0xFEEBE8E8,
        0x003E6F4B, 0x017AC0F9, 0xFE8B109C, 0x00E95AD6, 0xFE20125E, 0x00257688, 0x0002DE8A,
        0xFFC1FA76, 0xFF45BB33, 0x002BD31A, 0x01F875C6, 0x0389AF98, 0x0000527D,
    };

    bigint_word_t vc_actual[20];

    bigint_sub(vc_actual, va, vb, 20);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_sub_negative_result_20(void) {
    const bigint_word_t va[] = {
        0x01C17BF0, 0x01EC012D, 0x03992F9E, 0x01F11716, 0x02F0F7A3, 0x00FCEBFB, 0x03CEAE1A,
        0x02DB42EB, 0x021AE208, 0x01BA7181, 0x0042367E, 0x028BABFA, 0x0381C241, 0x0021E884,
        0x01D70AA5, 0x0352A584, 0x03423231, 0x00126B61, 0x006099A2, 0x0000073A,
    };
    const bigint_word_t vb[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vc_expected[] = {
        0x003FD868, 0x016AD258, 0x008472A6, 0xFFC2478E, 0x028DF5E2, 0xFFE8661B, 0x01141718,
        0xFFC190B5, 0xFE853F07, 0x0174EF64, 0xFF16A52A, 0x01DFEDA2, 0xFFDA8978, 0xFFFD2176,
        0x003E058A, 0x00BA44CD, 0xFFD42CE6, 0xFE078A3A, 0xFC765068, 0xFFFFAD83,
    };

    bigint_word_t vc_actual[20];

    bigint_sub(vc_actual, va, vb, 20);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_left_shift_255(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x03000000, 0x01D28831, 0x00A92F53, 0x03C8825C, 0x0014A168,
        0x00EB4708, 0x00ED0EFD, 0x0005831C, 0x00288DFA, 0x0228D38C, 0x00005C06,
    };

    bigint_word_t vb_actual[20];

    bigint_left_shift_255(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

#if 0
void test_bigint_right_shift_1(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x0128831C, 0x0292F537, 0x008825C2, 0x014A168F, 0x02B47080,
        0x02D0EFD3, 0x005831C3, 0x0288DFA0, 0x028D38C0, 0x0005C068,
    };

    bigint_word_t vb_actual[10];

    bigint_right_shift(vb_actual, va, 1);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_right_shift_17(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    const bigint_word_t vb_expected[] = {
        0x03D4DD28, 0x00970A92, 0x005A3C88, 0x01C2014A, 0x03BF4EB4,
        0x00C70ED0, 0x037E8058, 0x00E30288, 0x0301A28D, 0x00000005,
    };

    bigint_word_t vb_actual[10];

    bigint_right_shift(vb_actual, va, 17);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}
#endif

void test_bigint_right_shift_255(void) {
    const bigint_word_t va[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };

    const bigint_word_t vb_expected[] = {
        0x01722A82, 0x0177CB09, 0x00E71925, 0x0098E1DD, 0x0320A361,
        0x030C16EC, 0x01C0A974, 0x015C24FB, 0x01492750, 0x000B36FF,
    };

    bigint_word_t vb_actual[10];

    bigint_right_shift_255(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_right_shift_1(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x0128831C, 0x0292F537, 0x008825C2, 0x014A168F, 0x02B47080,
        0x02D0EFD3, 0x005831C3, 0x0288DFA0, 0x028D38C0, 0x0005C068,
    };

    bigint_word_t vb_actual[10];

    bigint_right_shift_1(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_by_div2_even(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x0128831C, 0x0292F537, 0x008825C2, 0x014A168F, 0x02B47080,
        0x02D0EFD3, 0x005831C3, 0x0288DFA0, 0x028D38C0, 0x0005C068,
    };

    bigint_word_t vb_actual[10];

    bigint_by_div2(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_by_div2_odd(void) {
    const bigint_word_t va[] = {
        0x02510639, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x01288313, 0x0292F537, 0x008825C2, 0x014A168F, 0x02B47080,
        0x02D0EFD3, 0x005831C3, 0x0288DFA0, 0x028D38C0, 0x0015C068,
    };

    bigint_word_t vb_actual[10];

    bigint_by_div2(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_mod_2_255(void) {
    const bigint_word_t va[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vb_expected[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1,
        0x011485E0, 0x02BA9702, 0x0319B236, 0x0395A301, 0x0005821D,
    };

    bigint_word_t vb_actual[10];

    bigint_mod_2_255(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_get_ith_bit_index_0_0(void) {
    bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    TEST_ASSERT_EQUAL(0, bigint_get_ith_bit(va, 0));
}

void test_bigint_get_ith_bit_index_0_1(void) {
    bigint_word_t va[] = {
        0x02510639, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    TEST_ASSERT_EQUAL(1, bigint_get_ith_bit(va, 0));
}

void test_bigint_get_ith_bit_index_254_0(void) {
    bigint_word_t va[] = {
        0x03FFFFED, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x000FFFFF,
    };

    TEST_ASSERT_EQUAL(0, bigint_get_ith_bit(va, 254));
}

void test_bigint_get_ith_bit_index_254_1(void) {
    bigint_word_t va[] = {
        0x03FFFFED, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x001FFFFF,
    };

    TEST_ASSERT_EQUAL(1, bigint_get_ith_bit(va, 254));
}

void test_bigint_csel_mask_0(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };

    bigint_word_t vc_actual[10];

    bigint_csel(vc_actual, va, vb, 0, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_csel_mask_FFFFFFFF(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vc_actual[10];

    bigint_csel(vc_actual, va, vb, 0xFFFFFFFF, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cadd_propagate_mask_0(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0288DF40, 0x028D38C0, 0x0005C068,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0042DAD9, 0x01F8A188, 0x000F9952,
    };
    const bigint_word_t vc_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0288DF40, 0x028D38C0, 0x0005C068,
    };

    bigint_word_t vc_actual[10];

    bigint_cadd_propagate(vc_actual, va, vb, 0, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cadd_propagate_mask_FFFFFFFF(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0288DF40, 0x028D38C0, 0x0005C068,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0042DAD9, 0x01F8A188, 0x000F9952,
    };
    const bigint_word_t vc_expected[] = {
        0x020C3F3F, 0x025B5967, 0x011064D1, 0x02F0061A, 0x0029069E,
        0x00B4F195, 0x02E1A2CC, 0x02CBBA19, 0x0085DA48, 0x001559BB,
    };

    bigint_word_t vc_actual[10];

    bigint_cadd_propagate(vc_actual, va, vb, 0xFFFFFFFF, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cgt_sub_much_greater(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x02C1376C, 0x03861208, 0x02FCCDDC, 0x015B7B89, 0x002C9F32,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x038FCECC, 0x019FD865, 0x02137DA8, 0x0138B194, 0x013C41CE,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vc_actual[10];

    TEST_ASSERT_EQUAL(1, bigint_cgt_sub(vc_actual, va, vb, 10));

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cgt_sub_much_smaller(void) {
    const bigint_word_t va[] = {
        0x02C1376C, 0x03861208, 0x02FCCDDC, 0x015B7B89, 0x002C9F32,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vc_expected[] = {
        0x02C1376C, 0x03861208, 0x02FCCDDC, 0x015B7B89, 0x002C9F32,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };

    bigint_word_t vc_actual[10];

    TEST_ASSERT_EQUAL(0, bigint_cgt_sub(vc_actual, va, vb, 10));

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cgt_sub_one_greater(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x02510637, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vc_expected[] = {
        0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };

    bigint_word_t vc_actual[10];

    TEST_ASSERT_EQUAL(1, bigint_cgt_sub(vc_actual, va, vb, 10));

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cgt_sub_one_smaller(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x02510639, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vc_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vc_actual[10];

    TEST_ASSERT_EQUAL(0, bigint_cgt_sub(vc_actual, va, vb, 10));

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_cgt_sub_equal(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vc_expected[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };

    bigint_word_t vc_actual[10];

    TEST_ASSERT_EQUAL(1, bigint_cgt_sub(vc_actual, va, vb, 10));

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_carry_propagate_positive(void) {
    bigint_word_t va[] = {
        0x12510638, 0x1125EA6A, 0x11104B81, 0x12942D1A, 0x1168E0FC,
        0x11A1DFA3, 0x10B06383, 0x1111BF3C, 0x111A717D, 0x000B80CD,
    };
    const bigint_word_t vb_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_carry_propagate(va, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, va, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_carry_propagate_negative(void) {
    bigint_word_t va[] = {
        0x11AEF9C8, 0x12DA158D, 0x12EFB476, 0x116BD2DD, 0x12971EFB,
        0x125E2054, 0x134F9C74, 0x12EE40BB, 0x02E58E7A, 0xFFF47F2E,
    };
    const bigint_word_t vb_expected[] = {
        0x01AEF9C8, 0x02DA1591, 0x02EFB47A, 0x016BD2E1, 0x02971EFF,
        0x025E2058, 0x034F9C78, 0x02EE40BF, 0x02E58E7E, 0xFFF47F2E,
    };

    bigint_carry_propagate(va, 10);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, va, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_carry_propagate_positive_double_size(void) {
    bigint_word_t va[] = {
        0x12510638, 0x1125EA6A, 0x11104B81, 0x12942D1A, 0x1168E0FC, 0x11A1DFA3, 0x10B06383,
        0x12237B3C, 0x1234E2FE, 0x11D7019E, 0x123BB38C, 0x131356EB, 0x13000190, 0x1385BD8B,
        0x136C0255, 0x1131311A, 0x126313F0, 0x12042DA9, 0x109F8A14, 0x0000F991,
    };
    const bigint_word_t vb_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100, 0x01A1DFA7, 0x00B06387,
        0x02237B40, 0x0234E302, 0x01D701A2, 0x023BB390, 0x031356EF, 0x03000194, 0x0385BD8F,
        0x036C0259, 0x0131311E, 0x026313F4, 0x02042DAD, 0x009F8A18, 0x0000F995,
    };

    bigint_carry_propagate(va, 20);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, va, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_mul(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };

    bigint_word_t vc_actual[20];

    bigint_mul(vc_actual, va, vb);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_mul_small(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb = 0x00000012;
    const bigint_word_t vc_expected[] = {
        0x01B26FF0, 0x00AA7BC6, 0x03254F5F, 0x026B2C20, 0x015FD20B,
        0x0161B9C4, 0x0066FF85, 0x033F7283, 0x03DBFB16, 0x00CF0EB6,
    };

    bigint_word_t vc_actual[10];

    bigint_mul_small(vc_actual, va, vb);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_barrett_Q_1(void) {
    const bigint_word_t vA[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vQ_expected[] = {
        0x01722A88, 0x0177CB09, 0x00E71925, 0x0098E1DD, 0x0320A361,
        0x030C16EC, 0x01C0A974, 0x015C24FB, 0x01492750, 0x000B36FF,
    };

    bigint_word_t vQ_actual[10], vR_actual[10];

    bigint_barrett(vQ_actual, vR_actual, vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vQ_expected, vQ_actual, sizeof(vQ_expected) / sizeof(vQ_expected[0]));
}

void test_bigint_barrett_R_1(void) {
    const bigint_word_t vA[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vR_expected[] = {
        0x00FACBA0, 0x00654087, 0x003B9ABE, 0x018792F4, 0x03CF21F7,
        0x02FA3972, 0x00072AAC, 0x00F070E0, 0x00038DF8, 0x001A9711,
    };

    bigint_word_t vQ_actual[10], vR_actual[10];

    bigint_barrett(vQ_actual, vR_actual, vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vR_actual, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_barrett_Q_2(void) {
    const bigint_word_t vA[] = {
        0x03FFFFEC, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x031FFFFF, 0x01D28831, 0x00A92F53, 0x03C8825C, 0x0014A168,
        0x00EB4708, 0x00ED0EFD, 0x0005831C, 0x00288DFA, 0x0228D38C, 0x00005C06,
    };
    const bigint_word_t vQ_expected[] = {
        0x0251063F, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vQ_actual[10], vR_actual[10];

    bigint_barrett(vQ_actual, vR_actual, vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vQ_expected, vQ_actual, sizeof(vQ_expected) / sizeof(vQ_expected[0]));
}

void test_bigint_barrett_R_2(void) {
    const bigint_word_t vA[] = {
        0x03FFFFEC, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x031FFFFF, 0x01D28831, 0x00A92F53, 0x03C8825C, 0x0014A168,
        0x00EB4708, 0x00ED0EFD, 0x0005831C, 0x00288DFA, 0x0228D38C, 0x00005C06,
    };
    const bigint_word_t vR_expected[] = {
        0x00037699, 0x01D06635, 0x00359AE4, 0x00FF593F, 0x02C8B30C,
        0x0303996B, 0x0117630C, 0x005131C3, 0x00F66C98, 0x001A8F88,
    };

    bigint_word_t vQ_actual[10], vR_actual[10];

    bigint_barrett(vQ_actual, vR_actual, vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vR_actual, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_reduce_25519_1(void) {
    bigint_word_t vA[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vR_expected[] = {
        0x00FACBA0, 0x00654087, 0x003B9ABE, 0x018792F4, 0x03CF21F7,
        0x02FA3972, 0x00072AAC, 0x00F070E0, 0x00038DF8, 0x001A9711,
    };

    bigint_reduce_25519(vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vA, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_reduce_25519_2(void) {
    bigint_word_t vA[] = {
        0x03FFFFEC, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x031FFFFF, 0x01D28831, 0x00A92F53, 0x03C8825C, 0x0014A168,
        0x00EB4708, 0x00ED0EFD, 0x0005831C, 0x00288DFA, 0x0228D38C, 0x00005C06,
    };
    const bigint_word_t vR_expected[] = {
        0x00037699, 0x01D06635, 0x00359AE4, 0x00FF593F, 0x02C8B30C,
        0x0303996B, 0x0117630C, 0x005131C3, 0x00F66C98, 0x001A8F88,
    };

    bigint_reduce_25519(vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vA, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_reduce_barrett_1(void) {
    bigint_word_t vA[] = {
        0x0181A388, 0x00812ED5, 0x0314BCF8, 0x022ECF88, 0x006301C1, 0x011485E0, 0x02BA9702,
        0x0319B236, 0x0395A301, 0x0045821D, 0x012B9154, 0x00ABBE58, 0x03A738C9, 0x0024C70E,
        0x0199051B, 0x029860B7, 0x036E054B, 0x020AE127, 0x03EA493A, 0x000059B7,
    };
    const bigint_word_t vR_expected[] = {
        0x00FACBA0, 0x00654087, 0x003B9ABE, 0x018792F4, 0x03CF21F7,
        0x02FA3972, 0x00072AAC, 0x00F070E0, 0x00038DF8, 0x001A9711,
    };

    bigint_reduce_barrett(vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vA, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_reduce_barrett_2(void) {
    bigint_word_t vA[] = {
        0x03FFFFEC, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF, 0x03FFFFFF,
        0x03FFFFFF, 0x03FFFFFF, 0x031FFFFF, 0x01D28831, 0x00A92F53, 0x03C8825C, 0x0014A168,
        0x00EB4708, 0x00ED0EFD, 0x0005831C, 0x00288DFA, 0x0228D38C, 0x00005C06,
    };
    const bigint_word_t vR_expected[] = {
        0x00037699, 0x01D06635, 0x00359AE4, 0x00FF593F, 0x02C8B30C,
        0x0303996B, 0x0117630C, 0x005131C3, 0x00F66C98, 0x001A8F88,
    };

    bigint_reduce_barrett(vA);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vR_expected, vA, sizeof(vR_expected) / sizeof(vR_expected[0]));
}

void test_bigint_to_montgomery(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t va_bar_expected[] = {
        0x0003769A, 0x01D06635, 0x00359AE4, 0x00FF593F, 0x02C8B30C,
        0x0303996B, 0x0117630C, 0x005131C3, 0x00F66C98, 0x001A8F88,
    };

    bigint_word_t va_bar_actual[10];

    bigint_to_montgomery(va_bar_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(va_bar_expected, va_bar_actual, sizeof(va_bar_expected) / sizeof(va_bar_expected[0]));
}

void test_bigint_montgomery_redc(void) {
    const bigint_word_t vc_bar_R[] = {
        0x02848DBE, 0x00AC83E8, 0x00B172AF, 0x0241C642, 0x017AF3E1, 0x00196FF2, 0x030A9516,
        0x030DDEA4, 0x02855B26, 0x0336BC19, 0x03E301F7, 0x00DBD5AC, 0x01F0926D, 0x0385DCD2,
        0x036AC865, 0x003A6C73, 0x02B8DE91, 0x014D79CF, 0x03CB2783, 0x00006F47,
    };
    const bigint_word_t vc_bar_expected[] = {
        0x029D1DFD, 0x0383CA09, 0x006C7C1B, 0x010FE81D, 0x005F855C,
        0x00924388, 0x00882AD2, 0x01D860A0, 0x0043896C, 0x00193643,
    };

    bigint_word_t vc_bar_actual[10];

    bigint_montgomery_redc(vc_bar_actual, vc_bar_R);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_bar_expected, vc_bar_actual, sizeof(vc_bar_expected) / sizeof(vc_bar_expected[0]));
}

void test_bigint_from_montgomery(void) {
    const bigint_word_t va_bar[] = {
        0x0003769A, 0x01D06635, 0x00359AE4, 0x00FF593F, 0x02C8B30C,
        0x0303996B, 0x0117630C, 0x005131C3, 0x00F66C98, 0x001A8F88,
    };
    const bigint_word_t va_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t va_actual[10];

    bigint_from_montgomery(va_actual, va_bar);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(va_expected, va_actual, sizeof(va_expected) / sizeof(va_expected[0]));
}

void test_bigint_modexp_25519_e_0(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_25519_e_1(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_25519_e_2(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000002, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x00DA6A3A, 0x00979F5A, 0x03B87A9A, 0x02F1CEA0, 0x01F04D44,
        0x037A6C7A, 0x02743719, 0x03A61177, 0x016304A0, 0x0002BD32,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_25519_e_3(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000003, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x02B30472, 0x002429EA, 0x0025002F, 0x00BDD99D, 0x00D2963E,
        0x002A7D0B, 0x032AC4CC, 0x02F213F2, 0x02613262, 0x001B7885,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_25519_e_2_254(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00100000,
    };
    const bigint_word_t vc_expected[] = {
        0x0214A916, 0x02516CF9, 0x00116D01, 0x030A95D9, 0x032BA0EF,
        0x00FE0A21, 0x00EC1BEB, 0x00DF0A91, 0x0353ED9E, 0x000A0A56,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_25519(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x0110D4F7, 0x03403B58, 0x013BFB7A, 0x024E82E2, 0x03B3D702,
        0x00160EFB, 0x01828B3B, 0x035EDD18, 0x03AE9A97, 0x0012B930,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_25519(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery_e_0(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery_e_1(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery_e_2(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000002, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x00DA6A3A, 0x00979F5A, 0x03B87A9A, 0x02F1CEA0, 0x01F04D44,
        0x037A6C7A, 0x02743719, 0x03A61177, 0x016304A0, 0x0002BD32,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery_e_3(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000003, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    };
    const bigint_word_t vc_expected[] = {
        0x02B30472, 0x002429EA, 0x0025002F, 0x00BDD99D, 0x00D2963E,
        0x002A7D0B, 0x032AC4CC, 0x02F213F2, 0x02613262, 0x001B7885,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery_e_2_254(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00100000,
    };
    const bigint_word_t vc_expected[] = {
        0x0214A916, 0x02516CF9, 0x00116D01, 0x030A95D9, 0x032BA0EF,
        0x00FE0A21, 0x00EC1BEB, 0x00DF0A91, 0x0353ED9E, 0x000A0A56,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_modexp_montgomery(void) {
    const bigint_word_t vb[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t ve[] = {
        0x03BB3907, 0x01356EF8, 0x0000194C, 0x005BD8FC, 0x02C0259E,
        0x031311ED, 0x02313F44, 0x0085B6D9, 0x03F14310, 0x001F32A4,
    };
    const bigint_word_t vc_expected[] = {
        0x0110D4F7, 0x03403B58, 0x013BFB7A, 0x024E82E2, 0x03B3D702,
        0x00160EFB, 0x01828B3B, 0x035EDD18, 0x03AE9A97, 0x0012B930,
    };

    bigint_word_t vc_actual[10];

    bigint_modexp_montgomery(vc_actual, vb, ve);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vc_expected, vc_actual, sizeof(vc_expected) / sizeof(vc_expected[0]));
}

void test_bigint_inv_exp(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x032C0DF7, 0x03241C65, 0x03120B48, 0x0220AB67, 0x01F16386,
        0x00EFBC26, 0x01C2BADC, 0x00128057, 0x01CB73A4, 0x0008D8EA,
    };

    bigint_word_t vb_actual[10];

    bigint_inv_exp(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

void test_bigint_inv_by(void) {
    const bigint_word_t va[] = {
        0x02510638, 0x0125EA6E, 0x01104B85, 0x02942D1E, 0x0168E100,
        0x01A1DFA7, 0x00B06387, 0x0111BF40, 0x011A7181, 0x000B80D1,
    };
    const bigint_word_t vb_expected[] = {
        0x032C0DF7, 0x03241C65, 0x03120B48, 0x0220AB67, 0x01F16386,
        0x00EFBC26, 0x01C2BADC, 0x00128057, 0x01CB73A4, 0x0008D8EA,
    };

    bigint_word_t vb_actual[10];

    bigint_inv_by(vb_actual, va);

    TEST_ASSERT_EQUAL_HEX32_ARRAY(vb_expected, vb_actual, sizeof(vb_expected) / sizeof(vb_expected[0]));
}

int main(int argc, char **argv) {
    utils_init();

    HAL_Delay(200);

    UNITY_BEGIN();
    RUN_TEST(test_bigint_add_10);
    RUN_TEST(test_bigint_add_20);
    RUN_TEST(test_bigint_sub_positive_result_10);
    RUN_TEST(test_bigint_sub_negative_result_10);
    RUN_TEST(test_bigint_sub_positive_result_20);
    RUN_TEST(test_bigint_sub_negative_result_20);
    RUN_TEST(test_bigint_left_shift_255);
    // RUN_TEST(test_bigint_right_shift_1);
    // RUN_TEST(test_bigint_right_shift_17);
    RUN_TEST(test_bigint_right_shift_255);
    RUN_TEST(test_bigint_right_shift_1);
    RUN_TEST(test_bigint_by_div2_even);
    RUN_TEST(test_bigint_by_div2_odd);
    RUN_TEST(test_bigint_mod_2_255);
    RUN_TEST(test_bigint_get_ith_bit_index_0_0);
    RUN_TEST(test_bigint_get_ith_bit_index_0_1);
    RUN_TEST(test_bigint_get_ith_bit_index_254_0);
    RUN_TEST(test_bigint_get_ith_bit_index_254_1);
    RUN_TEST(test_bigint_csel_mask_0);
    RUN_TEST(test_bigint_csel_mask_FFFFFFFF);
    RUN_TEST(test_bigint_cadd_propagate_mask_0);
    RUN_TEST(test_bigint_cadd_propagate_mask_FFFFFFFF);
    RUN_TEST(test_bigint_cgt_sub_much_greater);
    RUN_TEST(test_bigint_cgt_sub_much_smaller);
    RUN_TEST(test_bigint_cgt_sub_one_greater);
    RUN_TEST(test_bigint_cgt_sub_one_smaller);
    RUN_TEST(test_bigint_cgt_sub_equal);
    RUN_TEST(test_bigint_carry_propagate_positive);
    RUN_TEST(test_bigint_carry_propagate_negative);
    RUN_TEST(test_bigint_mul);
    RUN_TEST(test_bigint_mul_small);
    RUN_TEST(test_bigint_barrett_Q_1);
    RUN_TEST(test_bigint_barrett_R_1);
    RUN_TEST(test_bigint_barrett_Q_2);
    RUN_TEST(test_bigint_barrett_R_2);
    RUN_TEST(test_bigint_reduce_25519_1);
    RUN_TEST(test_bigint_reduce_25519_2);
    RUN_TEST(test_bigint_reduce_barrett_1);
    RUN_TEST(test_bigint_reduce_barrett_2);
    RUN_TEST(test_bigint_to_montgomery);
    RUN_TEST(test_bigint_montgomery_redc);
    RUN_TEST(test_bigint_from_montgomery);
    RUN_TEST(test_bigint_modexp_25519_e_0);
    RUN_TEST(test_bigint_modexp_25519_e_1);
    RUN_TEST(test_bigint_modexp_25519_e_2);
    RUN_TEST(test_bigint_modexp_25519_e_3);
    RUN_TEST(test_bigint_modexp_25519_e_2_254);
    RUN_TEST(test_bigint_modexp_25519);
    RUN_TEST(test_bigint_modexp_montgomery_e_0);
    RUN_TEST(test_bigint_modexp_montgomery_e_1);
    RUN_TEST(test_bigint_modexp_montgomery_e_2);
    RUN_TEST(test_bigint_modexp_montgomery_e_3);
    RUN_TEST(test_bigint_modexp_montgomery_e_2_254);
    RUN_TEST(test_bigint_modexp_montgomery);
    RUN_TEST(test_bigint_inv_exp);
    // RUN_TEST(test_bigint_inv_by);
    UNITY_END();

    while (1) {
    }
}
