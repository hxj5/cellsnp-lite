/* SNP operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_SNP_H
#define CSP_SNP_H

#include "htslib/sam.h"
#include "kvec.h"

/* 
* SNP List API
*/
/*@abstruct    The basic data structure to store SNP info.
@param chr     Name of chr.
@param pos     0-based coordinate in the reference sequence.
@param ref     Ref base (a letter). 0 means no ref in the input SNP file for the pos.
@param alt     Alt base (a letter). 0 means no alt in the input SNP file for the pos.
 */
typedef struct {
    char *chr;   
    hts_pos_t pos; 
    int8_t ref, alt;
} snp_t;

/*@abstract  Initilize the snp_t structure.
@return      Pointer to the snp_t structure if success, NULL if failure.
@note        The pointer returned successfully by snp_init() should be freed
             by snp_destroy() when no longer used.
 */
snp_t* snp_init(void);
void snp_destroy(snp_t *p); 
void snp_reset(snp_t *p); 

typedef kvec_t(snp_t*) snplist_t;  
#define snplist_destroy(v) {				\
    size_t __j;						\
    for (__j = 0; __j < kv_size(v); __j++) snp_destroy(kv_A(v, __j));	\
    kv_destroy(v);					\
    (v).a = NULL; (v).m = (v).n = 0;                    \
}

/*@abstract    Extract SNP info from bcf/vcf file.
@param fn      Filename of bcf/vcf.
@param pl      Pointer to client data used to store the extracted SNP info.
@param ret     Pointer to store running state. 0 if success, -1 otherwise.
@param print_skip  If print the skipped SNPs. 1, yes, 0, no.
@return        Num of elements successfully added into the snplist.

@note          If length of Ref or Alt is larger than 1, then the SNP would be skipped.
               If length of Ref or Alt is 0, then their values would be infered during pileup.
 */
size_t get_snplist_from_vcf(const char *fn, snplist_t *pl, int *ret, int print_skip);
#define get_snplist(fn, pl, ret, print_skip) get_snplist_from_vcf(fn, pl, ret, print_skip)

/*
 * Bi-Allele API
 */
typedef struct {
    int8_t ref, alt;
} biallele_t;

biallele_t* biallele_init(void);
void biallele_destroy(biallele_t *p);
void biallele_reset(biallele_t *p);

#endif

