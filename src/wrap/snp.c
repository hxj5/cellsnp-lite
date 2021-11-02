/* SNP operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "kvec.h"
#include "jstring.h"
#include "snp.h"

/* 
* SNP List API
*/

snp_t* snp_init(void) { return (snp_t*) calloc(1, sizeof(snp_t)); }

void snp_destroy(snp_t *p) { 
    if (p) { free(p); } 
}

void snp_reset(snp_t *p) {
    if (p) { memset(p, 0, sizeof(snp_t)); }
}

snp1_t* snp1_init(void) { return (snp1_t*) calloc(1, sizeof(snp1_t)); }

void snp1_destroy(snp1_t *p) { 
    if (p) { free(p); } 
}

chrsnp_t* chrsnp_init(void) { return (chrsnp_t*) calloc(1, sizeof(chrsnp_t)); }

void chrsnp_destroy(chrsnp_t *p) {
    if (p) { 
        if (p->chrom) { free(p->chrom); }
        if (p->is_raw) { snplist_destroy(p->sl); }
        free(p);
    }
}

void list_chrsnp_destroy(list_chrsnp_t p) {
    size_t i;
    for (i = 0; i < kv_size(p); i++) { chrsnp_destroy(kv_A(p, i)); }
    kv_destroy(p);
}

msnplist_t* msnplist_init(void) { return (msnplist_t*) calloc(1, sizeof(msnplist_t)); }
void msnplist_destroy(msnplist_t *p) {
    if (p) {
        list_chrsnp_destroy(p->cl);
        if (p->hcs) { kh_destroy(p->hcs); }
        free(p);
    }
}

size_t msnplist_size(msnplist_t *p) { return p->nsnp; }
size_t msnplist_nchrom(msnplist_t *p) { return p->nch; }

int msnplist_push(msnplist_t *pl, char *chrom, hts_pos_t pos, int8_t ref, int8_t alt) {
    chrsnp_t *cs = NULL;
    snp_t *ip = NULL;
    khiter_t k;
    int ret = -1, r, new_cs = 0;
    k = kh_get(chrsnp, pl->hcs, chrom);
    if (k == kh_end(pl->hcs)) {
        if (NULL == (cs = chrsnp_init())) { goto fail; }
        new_cs = 1;
        cs->chrom = safe_strdup(chrom);
        if (NULL == cs->chrom) { goto fail; }
        k = kh_put(chrsnp, pl->hcs, chrom, &r);
        kh_val(pl->hcs, k) = cs;
        kv_push(chrsnp_t*, pl->cl, cs);
        pl->nch++;
    } else { cs = kh_val(pl->hcs, k); }
    if (NULL == (ip = snp_init())) { goto fail; }
    ip->pos = pos; ip->ref = ref; ip->alt = alt;
    kv_push(snp_t*, cs->sl, ip);
    pl->nsnp++;
    return 0;
  fail:
    if (ip) { snp_destroy(ip); }
    if (cs && new_cs) { chrsnp_destroy(cs); }
    return ret;
}

static int msnplist_shrink(msnplist_t *pl) {
    chrsnp_t *cs;
    int i;
    for (i = 0; i < kv_size(pl->cl); i++) {
        cs = kv_A(pl->cl, i);
        kv_resize(snp_t*, cs->sl, kv_size(cs->sl));
    }
    return 0;
}

void msnplist_push_done(msnplist_t *pl) {
    size_t i;
    chrsnp_t *cs;
    msnplist_shrink(pl);
    for (i = 0; i < kv_size(pl->cl); i++) {
        cs = kv_A(pl->cl, i);
        if (cs->ibeg < 0) { cs->ibeg = 0; }
        if (cs->iend > kv_size(cs->sl)) { cs->iend = kv_size(cs->sl); }
        if (cs->ibeg >= cs->iend) {
            cs->ibeg = 0;
            cs->iend = kv_size(cs->sl);
        }
        cs->is_raw = 1;
    }
}

int msnplist_itr_next(msnplist_t *pl, snp1_t *s) {
    chrsnp_t *cs;
    snp_t *s0;
    while (1) {
        if (pl->ich >= kv_size(pl->cl)) { return 0; }
        cs = kv_A(pl->cl, pl->ich);
        if (pl->isnp >= cs->iend) { 
            pl->ich++;
            cs = kv_A(pl->cl, pl->ich);
            pl->isnp = cs->ibeg;
            continue;
        }
        if (pl->isnp < cs->ibeg) { pl->isnp = cs->ibeg; }
        s0 = kv_A(cs->sl, pl->isnp);
        pl->isnp++;
        break;
    }
    s->chrom = cs->chrom;
    s->data = s0;
    return 1;
}

void msnplist_itr_rewind(msnplist_t *pl) {
    chrsnp_t *cs;
    if (! pl || kv_size(pl->cl) <= 0) { return; }
    cs = kv_A(pl->cl, 0);
    pl->ich = 0;
    pl->isnp = cs->ibeg >= 0 ? cs->ibeg : 0;
}

msnplist_t** msnplist_split(msnplist_t *pl, unsigned n) {
    int i, c;
    size_t m, t, r, s, s1, beg;
    chrsnp_t *cs, *cs1 = NULL;
    msnplist_t **a, l = NULL;
    a = (msnplist_t**) calloc(n, sizeof(msnplist_t*));
    if (NULL == a) { goto fail; }
    if (n <= 0) { goto fail; }
    m = pl->nsnp / n;
    r = pl->nsnp - n * m;   // #remaining snp
    t = n < r ? m + 1 : m;  // #snp for each split
    for (i = 0, c = 0, s = s1 = 0, beg = 0; i < n; i++, s = s1 = 0) {
        if (NULL == (l = msnplist_init())) { goto fail; }
        while (s1 < t && c < pl->nch) {
            cs = kv_A(pl->cl, c);
            if (NULL == (cs1 = chrsnp_init())) { goto fail; }
            cs1->chrom = strdup(cs->chrom);
            cs1->sl = cs->sl;
            cs1->ibeg = beg; 
            cs1->iend = kv_size(cs->sl);
            //cs1->is_raw = 0;
            kv_push(chrsnp_t*, l->cl, cs1); cs1 = NULL;
            s += kv_size(cs->sl);
            s1 += kv_size(cs->sl) - beg;
            c++; beg = 0;
        }
        l->nch = kv_size(l->cl);
        l->nsnp = s;
        l->ich = 0;
        l->isnp = 0;
        if (s1 > t) {
            if (kv_size(l->cl) > 0) { 
                cs = kv_A(l->cl, kv_size(l->cl) - 1);
                cs->iend -= s1 - t;
                beg = cs->iend;
            } // else: should not come here!
        } else {
            c++;
            beg = 0;
        }
        a[i] = l; l = NULL;
    }
    return a;
  fail:
    int j;
    if (cs1) { chrsnp_destroy(cs1); }
    if (l) { msnplist_destroy(l); }
    if (a) {
        for (j = 0; j < i; j++) { 
            msnplist_destroy(a[j]);
        } free(a);
    }
    return NULL;
}

size_t get_snplist_from_vcf(const char *fn, msnplist_t *pl, int *ret, int print_skip) {
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    char *chrom = NULL;
    int8_t ref, alt;
    size_t l, m, n = 0;  
    int r;
    *ret = -1;
    if (NULL == fn || NULL == pl) { return 0; }
    if (NULL == (fp = hts_open(fn, "rb"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, fn); return 0; }
    if (NULL == (hdr = bcf_hdr_read(fp))) { fprintf(stderr, "[E::%s] could not read header for '%s'\n", __func__, fn); goto fail; }
    if (NULL == (rec = bcf_init1())) { fprintf(stderr, "[E::%s] could not initialize the bcf structure.\n", __func__); goto fail; }
    for (m = 1; (r = bcf_read1(fp, hdr, rec)) >= 0; m++) {
        chrom = bcf_hdr_id2name(hdr, rec->rid);
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele > 0) {
            if (1 == (l = strlen(rec->d.allele[0]))) { ref = rec->d.allele[0][0]; }
            else if (l > 1) { 
                if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: ref_len > 1.\n", __func__, m); }
                continue; 
            } // else: do nothing. keep ref = 0, 0 is its init value.               
            if (2 == rec->n_allele) {
                if (1 == (l = strlen(rec->d.allele[1]))) { alt = rec->d.allele[1][0]; }
                else if (l > 1) {
                    if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: alt_len > 1.\n", __func__, m); }
                    continue; 					
                } // else: do nothing. keep alt = 0, 0 is its init value.
            } else if (rec->n_allele > 2) {
                if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: n_allele > 2.\n", __func__, m); }
                continue;                 
            } // else: keep alt = 0.
        } // else: do nothing. keep ref = alt = 0, 0 is their init values.
        if (msnplist_push(pl, chrom, rec->pos, ref, alt) < 0) { goto fail; }
        n++;
    }
    if (-1 == r) { // end of bcf file.
        msnplist_push_done(pl);
    } else { 
        fprintf(stderr, "[E::%s] error when parsing '%s'\n", __func__, fn); 
        goto fail; 
    }
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    *ret = 0;
    return n;
  fail:
    if (rec) bcf_destroy1(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    return n;
}

/*
 * Bi-Allele API
 */
biallele_t* biallele_init(void) { return (biallele_t*) calloc(1, sizeof(biallele_t)); }
void biallele_destroy(biallele_t *p) { if (p) { free(p); } }
void biallele_reset(biallele_t *p) {}

