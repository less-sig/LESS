#!/usr/bin/env python3

q = 127
alignment = True
print("#include <stdint.h>")
print("const static uint32_t __fq127_lookup_table[] __attribute__((aligned(64))) = { ")

for i in range(q):
    l = [str((i * j) % q) for j in range(q)]
    print(", ".join(l), end='')
    if alignment:
        print(", 0,", end='')
    print()

print("};")
