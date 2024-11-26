#include <stdio.h>
#include "fq_arith.h"

int main() {
    printf("static const uint8_t inv_table[127] __attribute__((align(64))) = {\n0");
    for (uint32_t i = 1; i<127;i++) {
        printf(", %d", fq_inv(i));
    }

    printf("\n} ;\n");
}