/*******************************************************************************
 NIST-developed software is provided by NIST as a public service. You may use, 
 copy, and distribute copies of the software in any medium, provided that you 
 keep intact this entire notice. You may improve, modify, and create 
 derivative works of the software or any portion of the software, and you may 
 copy and distribute such modifications or works. Modified works should carry 
 a notice stating that you changed the software and should note the date and 
 nature of any such change. Please explicitly acknowledge the National 
 Institute of Standards and Technology as the source of the software.
 
 NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY 
 OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, 
 INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, 
 FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST 
 NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE 
 UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. 
 NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE 
 SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE 
 CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
  
 You are solely responsible for determining the appropriateness of using and 
 distributing the software and you assume all risks associated with its use, 
 including but not limited to the risks and costs of program errors, 
 compliance with applicable laws, damage to or loss of data, programs or 
 equipment, and the unavailability or interruption of operation. 
 This software is not intended to be used in any situation where a failure 
 could cause risk of injury or damage to property. 
 The software developed by NIST employees is not subject to copyright 
 protection within the United States.
*******************************************************************************/

#ifndef kat_nist_rng_h
#define kat_nist_rng_h


#include <string.h>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>

#define KAT_NIST_RNG_SUCCESS     ( 0)
#define KAT_NIST_RNG_BAD_MAXLEN  (-1)
#define KAT_NIST_RNG_BAD_OUTBUF  (-2)
#define KAT_NIST_RNG_BAD_REQ_LEN (-3)

typedef struct {
    unsigned char   buffer[16];
    unsigned int    buffer_pos;
    unsigned long   length_remaining;
    unsigned char   key[32];
    unsigned char   ctr[16];
} KAT_NIST_AES_XOF_struct;

typedef struct {
    unsigned char   Key[32];
    unsigned char   V[16];
    int             reseed_counter;
} KAT_NIST_AES256_CTR_DRBG_struct;


void
KAT_NIST_AES256_CTR_DRBG_Update(unsigned char *provided_data,
                                unsigned char *Key,
                                unsigned char *V);

int
KAT_NIST_seedexpander_init(KAT_NIST_AES_XOF_struct *ctx,
                           unsigned char *seed,
                           unsigned char *diversifier,
                           unsigned long maxlen);

int
KAT_NIST_seedexpander(KAT_NIST_AES_XOF_struct *ctx, 
                      unsigned char *x, unsigned long xlen);

void
KAT_NIST_randombytes_init(unsigned char *entropy_input,
                          unsigned char *personalization_string,
                          int security_strength);

int
KAT_NIST_randombytes(unsigned char *x, unsigned long long xlen);

#endif /* kat_nist_rng_h */
