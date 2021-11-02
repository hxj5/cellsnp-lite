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

typedef struct {
    hts_pos_t pos;    // 0-based coordinate in the reference sequence.
    int8_t ref, alt;  // Ref/Alt base. 0 means no ref/alt in the input SNP file for the pos.
} snp_t;

snp_t* snp_init(void);
void snp_destroy(snp_t *p); 
void snp_reset(snp_t *p); 

// Do not free the elements!
typedef struct {
    char *chrom;
    snp_t *data;
} snp1_t;

snp1_t* snp1_init(void);
void snp1_destroy(snp1_t *p); 

typedef kvec_t(snp_t*) snplist_t;  
#define snplist_destroy(v) {				\
    size_t __j;						\
    for (__j = 0; __j < kv_size(v); __j++) snp_destroy(kv_A(v, __j));	\
    kv_destroy(v);					\
    (v).a = NULL; (v).m = (v).n = 0;                    \
}

// use shared chrom string to save memory
typedef struct {
    char *chrom;
    snplist_t sl;
    size_t ibeg, iend;  // Iteration: index of begin (included) and end (excluded) elements.
    int is_raw;
} chrsnp_t;

chrsnp_t* chrsnp_init(void);
void chrsnp_destroy(chrsnp_t *p);

typedef kvec_t(chrsnp_t*) list_chrsnp_t;
void list_chrsnp_destroy(list_chrsnp_t p);

KHASH_MAP_INIT_INT(chrsnp, chrsnp_t*)
typedef khash_t(chrsnp) map_cs_t;

typedef struct {
    list_chrsnp_t cl;
    map_cs_t *hcs;    // For `push`
    size_t nch, nsnp; // total number of chrom and snp.
    size_t ich, isnp; // Iteration: index of next available chrom (ich) and snp in chrom (isnp).
    void *data;
} msnplist_t;

msnplist_t* msnplist_init(void);
void msnplist_destroy(msnplist_t *p);
size_t msnplist_size(msnplist_t *p);
size_t msnplist_nchrom(msnplist_t *p);

//@return 0 if success, -1 otherwise
int msnplist_push(msnplist_t *pl, char *chrom, hts_pos_t pos, int8_t ref, int8_t alt);
void msnplist_push_done(msnplist_t *pl);

//@abstract Iter next snp
//@return 1 if success, 0 otherwise
int msnplist_itr_next(msnplist_t *pl, snp1_t *s);

void msnplist_itr_rewind(msnplist_t *pl);

msnplist_t** msnplist_split(msnplist_t *pl, unsigned n);

/*!@func
@abstract   Extract SNP info from bcf/vcf file.
@param fn   Filename of bcf/vcf.
@param pl   Pointer to client data used to store the extracted SNP info.
@param ret  Pointer to store running state. 0 if success, -1 otherwise.
@param print_skip  If print the skipped SNPs. 1, yes, 0, no.
@return     Num of elements successfully added into the snplist.
@note       If length of Ref or Alt is larger than 1, then the SNP would be skipped.
            If length of Ref or Alt is 0, then their values would be infered during pileup.
 */
size_t get_snplist_from_vcf(const char *fn, msnplist_t *pl, int *ret, int print_skip);
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

