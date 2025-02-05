from code_utils import sample_monomial, KeyGen, build_full_generator_matrix, apply_monomial, inverse_perm, combine_perms, verify_rsp, compress
import random
import numpy as np
import galois
from CF_faster_impl import CF5_faster

k = 30
n = 60
q = 29
Fq = galois.GF(q); 

#sanity check = 1 if you want to do extra checks (e.g., permutation is correctly reverted)
sanity_checks = 0

A, A_prime, perm, coeffs = KeyGen(Fq, n, k)

print("sk:")
print("--> perm = ", perm)
print("--> coeffs = ",coeffs)
print(" ")
print("Non sys part of G:")
print(A)
print(" ")
print("Non sys part of G':")
print(A_prime)
print(" ")

#Sample another monomial, apply it to G and compute canonical form
ephemeral_perm, ephemeral_coeffs = sample_monomial(Fq, n)
full_G = build_full_generator_matrix(Fq, n, k, A)
ephemeral_A = apply_monomial(Fq, n, k, full_G, ephemeral_perm, ephemeral_coeffs)
cmt = CF5_faster(ephemeral_A, k, q, Fq)

print("cmt:")
print(cmt)
print(" ")

#Now, obtain same cmt from G_prime
perm_inv = inverse_perm(n, perm) #inverse of permutation in secret key
full_response_perm = combine_perms(n, perm_inv, ephemeral_perm) #this would be the full permutation for the response
rsp = compress(n, k, full_response_perm)
verify_cmt = verify_rsp(Fq, q, n, k, A_prime, rsp)


print("Is cmt vefified?",str(cmt) == str(verify_cmt))

#####SOME SANITY CHECKS, NO NEED TO LOOK AT THIS
if sanity_checks: #check permutation encoding is valid (we use only 1 coefficients for monomial)
    
    #apply only perm
    full_G = build_full_generator_matrix(Fq, n, k, A)
    perm_A = apply_monomial(Fq, n, k, full_G, perm, [Fq(1) for _ in range(n)])

    #now, apply inverse perm and verify result is the same
    full_perm_G = build_full_generator_matrix(Fq, n, k, perm_A)
    check_A = apply_monomial(Fq, n, k, full_perm_G, perm_inv, [Fq(1) for _ in range(n)])

    print("Is permutation encoding valid?",str(check_A) == str(A))
    
if sanity_checks: #verify combination of permutations is OK

    new_perm = combine_perms(n, perm_inv, ephemeral_perm)

    #apply only perm to G
    full_G = build_full_generator_matrix(Fq, n, k, A)
    new_A = apply_monomial(Fq, n, k, full_G, perm, [Fq(1) for _ in range(n)])

    #apply only perm
    ephemeral_A = apply_monomial(Fq, n, k, full_G, ephemeral_perm, [Fq(1) for _ in range(n)])

    #now, apply inverse perm and verify result is the same
    new_G = build_full_generator_matrix(Fq, n, k, new_A)
    check_A = apply_monomial(Fq, n, k, new_G, new_perm, [Fq(1) for _ in range(n)])

    print("Is combination of permutations valid?",str(check_A) == str(ephemeral_A))
