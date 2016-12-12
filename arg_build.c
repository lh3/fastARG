/* The MIT License

   Copyright (c) 2009, by Sanger Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <zlib.h>
#include "arg_data.h"
#include "kvec.h"
#include "khash.h"
#include "kbtree.h"

#define lt_alt(a, b) ((b) < (a))
KHASH_MAP_INIT_INT64(tract, int)
KHASH_SET_INIT_INT64(set64)

/***********************
 * defines and structs *
 ***********************/

#define AB_FHAP      0x01
#define AB_FSHOW_SID 0x02
#define AB_FSHOW_SEQ 0x04
#define AB_FBACKWARD 0x08
#define AB_FLONGEST  0x10

typedef kvec_t(int) intvec_t;

#define btcmp(a, b) (((a) < (b)) - ((b) < (a)))
KBTREE_INIT(bt, uint64_t, btcmp)

typedef struct {
	int beg, end, nid, z;
	int *s; // IMPORTANT: 0->hom1, 1->hom2, 2->N, >=3 for het (different from arg_data_t::seq)
	intvec_t *mut;
	kbtree_t(bt) *bt;
} arg_S1_t;

KHASH_MAP_INIT_INT(S, arg_S1_t)

typedef struct {
	int m, n, sid, nid;
	khash_t(S) *S;
	khash_t(tract) *tract;
	khash_t(set64) *C_set; // pairs where coalescence can be applied
	int has_M, *M_set; // sites where mutations can be applied
	int *count; // dimension: m*4
	int flag;
	int *phase[3], max_het;
	double ratio;
} arg_build_t;

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static int g_n_recombs = 0;

/*******************
 * Phasing related *
 *******************/

static int gen_phase_arrays(arg_build_t *ab, int **seq)
{
	int i, j, l, x, max_x;
	x = 3; max_x = 4;
	for (l = 0; l < 3; ++l)
		ab->phase[l] = (int*)calloc(max_x, sizeof(int));
	for (i = 0; i < ab->m; ++i) {
		for (j = 0; j < ab->n; ++j) {
			int c = seq[j][i]; // this is not cache efficient, but this is for future cache efficiency
			if (c < 3) continue; // phased or missing
			if (x == max_x) {
				max_x <<= 1;
				for (l = 0; l < 3; ++l)
					ab->phase[l] = (int*)realloc(ab->phase[l], max_x * sizeof(int));
			}
			seq[j][i] = c = x++;
			if (j%2 == 1) {
				ab->phase[0][c] = c-1;
				ab->phase[0][c-1] = c;
			}
			ab->phase[1][c] = c;
			ab->phase[2][c] = j;
		}
	}
	return x;
}

static int unphased_unphased(arg_build_t *ab, int sid, int c1, int c2, khash_t(set64) *set)
{
	int x, y, d1, d2, z;
	assert(c1 >= 3 && c2 >= 3);
	d1 = ab->phase[0][c1]; d2 = ab->phase[0][c2];

	// update phase[0] in the d2-loop
	if (d2) {
		x = d2;
		do {
			if (ab->phase[0][x] == c2) ab->phase[0][x] = c1;
			x = ab->phase[1][x];
		} while (x != d2);
	}

	// check if c2 is in the c1-loop
	x = c1; y = -1;
	do {
		if (x == c2) break;
		y = x;
		x = ab->phase[1][x];
	} while (x != c1);
	if (x == c2) { // in the c1-loop
		ab->phase[1][y] = ab->phase[1][c2]; // remove c2 from the c1-loop
	} else { // c2 is not in the c1-loop
		// merge d2-loop to d1-loop
		if (d1 && d2) {
			x = d1; do { y = x; x = ab->phase[1][x]; } while (x != d1);
			x = d2; do { z = x; x = ab->phase[1][x]; } while (x != d2);
			ab->phase[1][y] = d2; ab->phase[1][z] = d1;
		} else if (d1 == 0) d1 = ab->phase[0][c1] = d2;
		// merge c2-loop to c1-loop
		if (ab->phase[1][c2] != c2) {
			x = c1; do { y = x; x = ab->phase[1][x]; } while (x != c1);
			x = ab->phase[1][c2];
			while (x != c2) { // excluding c2
				z = x; x = ab->phase[1][x];
			}
			ab->phase[1][y] = ab->phase[1][c2]; ab->phase[1][z] = c1;
		}
	}

	// remove c2 and update phase[2]
	ab->phase[1][c2] = ab->phase[2][c2] = -1;
	ab->phase[2][c1] = sid;

	// put c1-loop and d1-loop to "set"; TODO: this can be improved
	if (d1) {
		x = d1;
		do {
			kh_put(set64, set, ab->phase[2][x], &z);
			x = ab->phase[1][x];
		} while (x != d1);
	}
	x = c1;
	do {
		kh_put(set64, set, ab->phase[2][x], &z);
		x = ab->phase[1][x];
	} while (x != c1);

	return c1;
}

static void phased_unphased(arg_build_t *ab, int i, int c1, int c2, khash_t(set64) *set)
{
	int x;
	assert(c1 <= 1 && c2 >= 3);

#define __phase(_c, _x) do {					\
		arg_S1_t *p;							\
		khint_t k;								\
		int ret, j = ab->phase[2][_x];			\
		k = kh_get(S, ab->S, j);				\
		assert(k != kh_end(ab->S));				\
		p = &kh_val(ab->S, k);					\
		p->s[i - p->beg] = (_c);				\
		++ab->count[i<<2|(_c)];					\
		--ab->count[i<<2|3];					\
		kh_put(set64, set, j, &ret);			\
		ab->phase[2][_x] = -1;					\
		j = _x;									\
		_x = ab->phase[1][_x];					\
		ab->phase[1][j] = -1;					\
	} while (0)

	// phase c2-loop
	x = c2; do __phase(c1, x); while (x != c2);
	if (ab->phase[0][c2]) { // then phase phase[0][c2]-loop
		x = ab->phase[0][c2];
		do __phase(1-c1, x); while (x != ab->phase[0][c2]);
	}
}

// requirement: c1>2 && c2>2
static inline int is_coalescable(const arg_build_t *ab, int c1, int c2)
{
	int d1 = ab->phase[0][c1], x;
	x = d1;
	do {
		if (x == c2) return 0;
		x = ab->phase[1][x];
	} while (x != d1);
	return 1;
}

/********************
 * find first tract *
 ********************/

static int find_tract(const arg_build_t *ab, const arg_S1_t *s1, const arg_S1_t *s2, int *is_coal)
{
	int i, beg, end, first;

	// check overlap; sister checking should be done by the caller
	*is_coal = 0;
	beg = s1->beg > s2->beg? s1->beg : s2->beg;
	end = s1->end < s2->end? s1->end : s2->end;
	if (beg >= end) return 0; // no overlap

#define __ft_loop2(_i) {							\
		int c1 = s1->s[(_i) - s1->beg];				\
		int c2 = s2->s[(_i) - s2->beg];				\
		if (c1 < 2 && c2 < 2) {						\
			if (c1 != c2) break;					\
		} else if (c1 > 2 && c2 > 2) {				\
			if (!is_coalescable(ab, c1, c2)) break;	\
		}											\
	}

	// get the first tract
	if (ab->flag & AB_FBACKWARD) {
		for (i = end - 1; i >= beg; --i) __ft_loop2(i);
		first = end - 1 - i;
	} else {
		for (i = beg; i < end; ++i) __ft_loop2(i);
		first = i - beg;
	}
	if (first == end - beg) *is_coal = 1;
	return first;
}

/********************************************
 * print sequences (for debugging purposes) *
 ********************************************/

static void print_seq1(const arg_build_t *ab, int sid, int beg, int end)
{
	int i;
	khint_t k;
	arg_S1_t *p;
	k = kh_get(S, ab->S, sid);
	p = &kh_val(ab->S, k);
	printf("Q\t%d\t%d\t[%d,%d]\t", kh_key(ab->S, k), p->nid, beg, end);
	for (i = beg; i < end; ++i) {
		if (i < p->beg || i >= p->end) putchar('.');
		else putchar(p->s[i - p->beg] > 3? '3' : '0'+p->s[i - p->beg]);
	}
	putchar('\n');
}

static void print_seq(const arg_build_t *ab)
{
	khint_t k;
	for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k)
		if (kh_exist(ab->S, k))
			print_seq1(ab, kh_key(ab->S, k), 0, ab->m);
}

/**************************************
 * Auxiliary functions for add/remove *
 **************************************/

static void update_M_set(arg_build_t *ab)
{
	int j;
	khint_t k;
	ab->has_M = 0;
	for (j = 0; j < ab->m; ++j) {
		int *c = ab->count + (j<<2);
		if (c[2] == 0 && c[3] == 0 && ((c[0] == 1 && c[1] >= 1) || (c[1] == 1 && c[0] >= 1))) { // NB: c[2] cannot be 1
			int x;
			if (c[0] == 1 && c[1] == 1) x = (drand48() < 0.5);
			else x = c[0] == 1? 0 : 1;
			for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k) {
				if (kh_exist(ab->S, k)) {
					arg_S1_t *p = &kh_val(ab->S, k);
					if (j >= p->beg && j < p->end && p->s[j - p->beg] == x) break;
				}
			}
			ab->M_set[j] = kh_key(ab->S, k);
			ab->has_M = 1;
		} else ab->M_set[j] = -1;
	}
}

// requirement: s corresponds to p; t to q
static inline void remove_from_bt(int first, int s, arg_S1_t *p, int t, arg_S1_t *q, int flag)
{
	uint64_t r;
	// remove from (s, p)
	r = (uint64_t)first<<32 | t;
	if (kb_get(bt, p->bt, r)) {
		kb_del(bt, p->bt, r);
		if (flag & AB_FLONGEST) __kb_get_first(uint64_t, p->bt, p->z);
		else p->z -= first;
	}
	// remove from (t, q)
	r = (uint64_t)first<<32 | s;
	if (kb_get(bt, q->bt, r)) {
		kb_del(bt, q->bt, r);
		if (flag & AB_FLONGEST) __kb_get_first(uint64_t, q->bt, q->z);
		else q->z -= first;
	}
}

// requirement: t absent from p->bt and s from q->bt (remove_from_bt
// should be called before hand if necessary)
static inline void insert_to_bt(int first, int s, arg_S1_t *p, int t, arg_S1_t *q, int flag)
{
	uint64_t r;
	// add to (s, p)
	r = (uint64_t)first<<32 | t;
	kb_put(bt, p->bt, r);
	if (flag & AB_FLONGEST) __kb_get_first(uint64_t, p->bt, p->z);
	p->z += first;
	// add to (t, q)
	r = (uint64_t)first<<32 | s;
	kb_put(bt, q->bt, r);
	if (flag & AB_FLONGEST) __kb_get_first(uint64_t, q->bt, q->z);
	q->z += first;
}

static void update_tract(arg_build_t *ab, int sid)
{
	khint_t iter, k;
	arg_S1_t *p;

	k = kh_get(S, ab->S, sid);
	assert(k != kh_end(ab->S));
	p = &kh_val(ab->S, k);

	for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k) {
		if (kh_exist(ab->S, k) && kh_key(ab->S, k) != sid) {
			int first_tract, old_tract;
			uint64_t key;
			arg_S1_t *q = &kh_val(ab->S, k);
			int is_coal = 0, ret;
			key = kh_key(ab->S, k) < sid? (uint64_t)kh_key(ab->S, k)<<32|sid : (uint64_t)sid<<32|kh_key(ab->S, k);
			first_tract = find_tract(ab, p, q, &is_coal);
			if (first_tract == 0) {
				iter = kh_get(tract, ab->tract, key);
				if (iter != kh_end(ab->tract)) {
					remove_from_bt(kh_val(ab->tract, iter), sid, p, kh_key(ab->S, k), q, ab->flag);
					// remove from ab->tract
					kh_del(tract, ab->tract, iter);
				}
				// remove from ab->C_set
				kh_del(set64, ab->C_set, kh_get(set64, ab->C_set, key));
				continue;
			}
		    if (is_coal) {
				kh_put(set64, ab->C_set, key, &ret);
				continue;
			}
			kh_del(set64, ab->C_set, kh_get(set64, ab->C_set, key)); // not coalescable
			iter = kh_put(tract, ab->tract, key, &ret);
			old_tract = kh_val(ab->tract, iter);
			kh_val(ab->tract, iter) = first_tract;
			if (ret == 0) remove_from_bt(old_tract, sid, p, kh_key(ab->S, k), q, ab->flag);
			insert_to_bt(first_tract, sid, p, kh_key(ab->S, k), q, ab->flag);
		}
	}
}

/************************
 * Add/remove sequences *
 ************************/

static int add_seq(arg_build_t *ab, const arg_S1_t *s0)
{
	int ret, sid, i;
	khint_t k;
	arg_S1_t *p;

	sid = ab->sid++;
	k = kh_put(S, ab->S, sid, &ret);
	p = &kh_val(ab->S, k);
	*p = *s0;
	for (i = p->beg; i < p->end; ++i) {
		int c = p->s[i - p->beg];
		if (c >= 3) {
			ab->phase[2][c] = sid;
			c = 3;
		}
		++ab->count[i<<2 | c];
	}
	update_M_set(ab); // update ->M_set
	update_tract(ab, sid); // update ->tract and ->C_set, O(mn)
	return sid;
}

static void remove_seq(arg_build_t *ab, int sid)
{
	khint_t k, iter;
	int j;
	arg_S1_t *p;

	// update ->S and ->count, ~O(1)
	k = kh_get(S, ab->S, sid);
	assert(k < kh_end(ab->S));
	p = &kh_value(ab->S, k);
	for (j = p->beg; j < p->end; ++j) {
		int c = p->s[j - p->beg];
		if (c >= 3) c = 3;
		--ab->count[j<<2 | c];
	}
	kh_del(S, ab->S, k); // after this, p is still a valid pointer pointing to the right place
	update_M_set(ab); // update ->M_set

	// update ->tract and ->C_set, O(m)
	for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k) {
		if (kh_exist(ab->S, k)) {
			uint64_t key = kh_key(ab->S, k) < sid? (uint64_t)kh_key(ab->S, k)<<32|sid : (uint64_t)sid<<32|kh_key(ab->S, k);
			iter = kh_get(tract, ab->tract, key);
			if (iter != kh_end(ab->tract)) {
				remove_from_bt(kh_val(ab->tract, iter), sid, p, kh_key(ab->S, k), &kh_val(ab->S, k), ab->flag);
				kh_del(tract, ab->tract, iter);
			}
			kh_del(set64, ab->C_set, kh_get(set64, ab->C_set, key));
		}
	}

	// destroy
	free(p->s);
	if (p->mut) free(p->mut->a);
	free(p->mut);
	kb_destroy(bt, p->bt); p->bt = 0;
}

/*********************************************
 * initialize and destroy arg_build_t struct *
 *********************************************/

static arg_build_t *init_build(const arg_data_t *ad, int flag, double ratio)
{
	arg_build_t *ab;
	int i, j, **abseq;

	// initialization
	ab = (arg_build_t*)calloc(1, sizeof(arg_build_t));
	ab->flag = flag;
	ab->ratio = ratio;
	if (ad->is_phased) ab->flag |= AB_FHAP;
	else ab->flag &= ~AB_FHAP;
	ab->m = ad->m; ab->sid = ab->nid = 0;
	ab->n = ad->n;
	ab->S = kh_init(S);
	ab->count = (int*)calloc(ab->m * 4, sizeof(int));
	ab->M_set = (int*)calloc(ab->m, sizeof(int));
	ab->tract = kh_init(tract);
	ab->C_set = kh_init(set64);
	for (i = 0; i < ab->m; ++i) ab->M_set[i] = -1;
	abseq = (int**)calloc(ab->n, sizeof(int*));

	// add sequences
	for (i = 0; i < ab->n; ++i) {
		char *seq = ad->seq[i];
		int *p;
		p = abseq[i] = (int*)calloc(ab->m, sizeof(int));
		for (j = 0; j < ad->m; ++j) {
			int c = seq[j];
			assert(c >= 0 && c <= 3);
			if (c == 3) c = 2; // missing data
			else if (c == 2) c = 3; // unphased data
			p[j] = c;
		}		
	}

	if (!(ab->flag & AB_FHAP)) ab->max_het = gen_phase_arrays(ab, abseq);
	for (i = 0; i < ab->n; ++i) {
		arg_S1_t s;
		s.beg = 0; s.end = ab->m; s.nid = ab->nid++; s.z = 0;
		s.mut = 0; s.bt = kb_init(bt, KB_DEFAULT_SIZE);
		s.s = abseq[i];
		add_seq(ab, &s);
	}
	printf("N\t%d\t%d\n", kh_size(ab->S), ab->m);
	free(abseq);
	return ab;
}

static void destroy_build(arg_build_t *ab)
{
	khint_t k;
	if (ab == 0) return;
	for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k) {
		if (kh_exist(ab->S, k)) {
			arg_S1_t *p = &kh_val(ab->S, k);
			kb_destroy(bt, p->bt);
			free(p->s);
			if (p->mut) kv_destroy(*p->mut);
		}
	}
	kh_destroy(S, ab->S);
	free(ab->count);
	free(ab->M_set);
	free(ab->phase[0]); free(ab->phase[1]); free(ab->phase[2]);
	kh_destroy(set64, ab->C_set);
	kh_destroy(tract, ab->tract);
	free(ab);
}

/*****************************
 * coalesce/mutate/recombine *
 *****************************/

static inline void print_mut(intvec_t *mut)
{
	if (mut == 0 || mut->n == 0) printf("0\n");
	else {
		int k;
		printf("%d", (int)mut->n);
		for (k = 0; k < mut->n; ++k)
			printf("\t%d", mut->a[k]);
		putchar('\n');
	}
}

static int coalesce(arg_build_t *ab, int s1, int s2)
{
	khint_t k;
	arg_S1_t s, *p[2];
	int i;
	uint64_t key;
	khash_t(set64) *set;

	if (ab->flag & AB_FSHOW_SID) printf("c\t(%d,%d) -> %d\n", s1, s2, ab->sid);
	key = s1 < s2? (uint64_t)s1<<32|s2 : (uint64_t)s2<<32|s1;
	assert(kh_get(set64, ab->C_set, key) != kh_end(ab->C_set));

	set = kh_init(set64);
	k = kh_get(S, ab->S, s1); p[0] = &kh_val(ab->S, k);
	k = kh_get(S, ab->S, s2); p[1] = &kh_val(ab->S, k);

	s.beg = p[0]->beg < p[1]->beg? p[0]->beg : p[1]->beg;
	s.end = p[0]->end > p[1]->end? p[0]->end : p[1]->end;
	s.nid = ab->nid++; s.z = 0;
	s.mut = 0; s.bt = kb_init(bt, KB_DEFAULT_SIZE);
	s.s = (int*)calloc(s.end - s.beg, sizeof(int));

	for (i = s.beg; i < s.end; ++i) {
		if (i >= p[0]->beg && i < p[0]->end) {
			if (i >= p[1]->beg && i < p[1]->end) {
				int c[2], x;
				c[0] = p[0]->s[i - p[0]->beg];
				c[1] = p[1]->s[i - p[1]->beg];
				x = c[0] < c[1]? 0 : 1;
				if (c[x] == 2) x = 1 - x;
				s.s[i - s.beg] = c[x];
				if (c[x] < 2 && c[1-x] >= 3) phased_unphased(ab, i, c[x], c[1-x], set);
				else if (c[x] >= 3 && c[1-x] >= 3)
					s.s[i - s.beg] = unphased_unphased(ab, ab->sid, c[x], c[1-x], set);
			} else s.s[i - s.beg] = p[0]->s[i - p[0]->beg]; // take p[0]->s
		} else s.s[i - s.beg] = p[1]->s[i - p[1]->beg]; // take p[1]->s
	}

	for (k = kh_begin(set); k < kh_end(set); ++k)
		if (kh_exist(set, k) && (int)kh_key(set, k) < ab->sid)
			update_tract(ab, (int)kh_key(set, k));
	kh_destroy(set64, set);

	printf("C\t%d\t%d\t%d\t%d\t", s.nid, p[0]->nid, p[0]->beg, p[0]->end); print_mut(p[0]->mut);
	printf("C\t%d\t%d\t%d\t%d\t", s.nid, p[1]->nid, p[1]->beg, p[1]->end); print_mut(p[1]->mut);
	remove_seq(ab, s1);
	remove_seq(ab, s2);
	return add_seq(ab, &s);
}

static int mutate(arg_build_t *ab, int s0)
{
	arg_S1_t s, *p;
	khint_t k;
	int i;

	if (ab->flag & AB_FSHOW_SID) printf("m\t%d -> %d\n", s0, ab->sid);
	k = kh_get(S, ab->S, s0); p = &kh_val(ab->S, k);
	s.beg = p->beg; s.end = p->end; s.nid = p->nid; s.z = 0;
	s.s = (int*)calloc(s.end - s.beg, sizeof(int));
	memcpy(s.s, p->s, (s.end - s.beg) * sizeof(int));
	s.mut = (intvec_t*)calloc(1, sizeof(intvec_t));
	s.bt = kb_init(bt, KB_DEFAULT_SIZE);
	if (p->mut)
		for (i = 0; i < p->mut->n; ++i)
			kv_push(int, s.mut[0], p->mut->a[i]);
	for (i = 0; i < ab->m; ++i) {
		if (ab->M_set[i] == s0) {
			assert(s.s[i - s.beg] < 2);
			s.s[i - s.beg] = 1 - s.s[i - s.beg];
			kv_push(int, s.mut[0], i);
		}
	}

	remove_seq(ab, s0);
	return add_seq(ab, &s);
}

static int recombine_aux(arg_build_t *ab, int s0)
{
	arg_S1_t s, *p;
	khint_t k;

	if (ab->flag & AB_FSHOW_SID) printf("r\t%d -> %d\n", s0, ab->sid);
	k = kh_get(S, ab->S, s0); p = &kh_val(ab->S, k);
	s.beg = p->beg; s.end = p->end; s.nid = ab->nid++; s.z = 0;
	s.s = (int*)calloc(s.end - s.beg, sizeof(int));
	memcpy(s.s, p->s, (s.end - s.beg) * sizeof(int));
	s.mut = 0; s.bt = kb_init(bt, KB_DEFAULT_SIZE);

	printf("R\t%d\t%d\t%d\t%d\t", s.nid, p->nid, p->beg, p->end); print_mut(p->mut);
	++g_n_recombs;
	remove_seq(ab, s0);
	return add_seq(ab, &s);
}

// recombination happens between x-1 and x
static int recombine(arg_build_t *ab, int x, int s0, int is_end)
{
	khint_t k;
	arg_S1_t z[2], *p;
	int r1, r2, s1;

	s1 = recombine_aux(ab, s0);
	if (ab->flag & AB_FSHOW_SID) printf("s\t%d -> (%d,%d) @ %d\n", s1, ab->sid, ab->sid+1, x);

	k = kh_get(S, ab->S, s1);
	p = &kh_val(ab->S, k);
	assert(x > p->beg && x < p->end);
	z[0].beg = p->beg; z[0].end = x; z[0].nid = p->nid;
	z[1].beg = x; z[1].end = p->end; z[1].nid = p->nid;
	z[0].z = z[1].z = 0;
	z[0].mut = z[1].mut = 0;
	z[0].bt = kb_init(bt, KB_DEFAULT_SIZE);
	z[1].bt = kb_init(bt, KB_DEFAULT_SIZE);

	z[0].s = (int*)calloc(z[0].end - z[0].beg, sizeof(int));
	z[1].s = (int*)calloc(z[1].end - z[1].beg, sizeof(int));
	memcpy(z[0].s, p->s, (x - p->beg) * sizeof(int));
	memcpy(z[1].s, p->s + (x - p->beg), (p->end - x) * sizeof(int));

	remove_seq(ab, s1);
	r1 = add_seq(ab, &z[0]);
	r2 = add_seq(ab, &z[1]);
	return is_end? r2 : r1;
}

/****************************************
 * pick a pair of sequence to recombine *
 ****************************************/

static int pick_pair(const arg_build_t *ab, int s[2])
{
	int ret = 0, sid = -1;
	khint_t k;
	arg_S1_t *p = 0;
	kbtree_t(bt) *bt;
	bt = kb_init(bt, KB_DEFAULT_SIZE);

	for (k = kh_begin(ab->S); k < kh_end(ab->S); ++k)
		if (kh_exist(ab->S, k))
			kb_put(bt, bt, (uint64_t)kh_val(ab->S, k).z<<32 | kh_key(ab->S, k));

#define bt_traverse_cal_sum(_r) {						\
		uint64_t r = *(_r);								\
		if (max_l == 0) max_l = r>>32;					\
		if (i >= thres && max_l > (int)(r>>32)) break;	\
		++i; sum += r>>32;								\
	}

#define bt_traverse_random(_r) {					\
		uint64_t r = *(_r);							\
		sum += r>>32;								\
		ret = r>>32; sid = (int)r;					\
		if ((double)sum/tot >= rand_num) break;		\
	}

#define bt_traverse_max(_r) {									\
		uint64_t r = *(_r);										\
		uint64_t cval;											\
		if (max == 0) max = r>>32;								\
		if (max != (int)(r>>32)) break;							\
		cval = r>>32<<32 | ((uint32_t)lrand48());				\
		if (cval > max_cval) {									\
			max_cval = cval; sid = (uint32_t)r; ret = r>>32;	\
		}														\
	}

	if (ab->ratio == 0.0) {
		int max = 0;
		uint64_t max_cval = 0;
		__kb_traverse(uint64_t, bt, bt_traverse_max);
	} else {
		int sum = 0, max_l = 0, i = 0, tot;
		double rand_num = drand48();
		int thres = (int)(kh_size(ab->S) * ab->ratio + .499);
		__kb_traverse(uint64_t, bt, bt_traverse_cal_sum);
		tot = sum; sum = 0;
		__kb_traverse(uint64_t, bt, bt_traverse_random);
	}
	assert(sid >= 0);
	s[0] = sid;
	k = kh_get(S, ab->S, sid);
	p = &kh_val(ab->S, k);
	sid = -1;
	if (ab->ratio == 0.0) {
		int max = 0;
		uint64_t max_cval = 0;
		__kb_traverse(uint64_t, p->bt, bt_traverse_max);
	} else {
		int sum = 0, max_l = 0, i = 0, tot;
		double rand_num = drand48();
		int thres = (int)(kb_size(p->bt) * ab->ratio + .499);
		__kb_traverse(uint64_t, p->bt, bt_traverse_cal_sum);
		tot = sum; sum = 0;
		__kb_traverse(uint64_t, p->bt, bt_traverse_random);
	}
	s[1] = sid;
	assert(sid >= 0);
	kb_destroy(bt, bt);
	return ret;
}

/*****************
 * the core loop *
 *****************/

void arg_build(arg_build_t *ab)
{
	khint_t k;
	int cnt = 0;
	while (kh_size(ab->S) > 1) {
		if (cnt % 1000 == 0)
			fprintf(stderr, "[arg_build] %d actions; %d sequences remain\n", cnt, (int)kh_size(ab->S));
		++cnt;
		if (ab->flag & AB_FSHOW_SEQ) print_seq(ab);
		if (kh_size(ab->C_set)) { // then coalesce
			for (k = kh_begin(ab->C_set); k < kh_end(ab->C_set); ++k)
				if (kh_exist(ab->C_set, k)) {
					coalesce(ab, kh_key(ab->C_set, k)>>32, (uint32_t)kh_key(ab->C_set, k));
					break;
				}
		} else if (ab->has_M) { // then mutate
			int i;
			for (i = 0; i < ab->m; ++i)
				if (ab->M_set[i] >= 0) {
					mutate(ab, ab->M_set[i]); break;
				}
		} else { // then recombine
			int s[2], beg, end, x, tract;
			arg_S1_t *p[2];
			tract = pick_pair(ab, s);
			assert(tract);
			k = kh_get(S, ab->S, s[0]); assert(k < kh_end(ab->S)); p[0] = &kh_val(ab->S, k);
			k = kh_get(S, ab->S, s[1]); assert(k < kh_end(ab->S)); p[1] = &kh_val(ab->S, k);
			if (ab->flag & AB_FBACKWARD) {
				end = p[0]->end < p[1]->end? p[0]->end : p[1]->end;
				beg = end - tract;
				assert(beg != 0);
				if (end == p[0]->end) {
					if (end == p[1]->end) x = (drand48() < 0.5);
					else x = 0;
				} else x = 1;
				coalesce(ab, recombine(ab, beg, s[x], 1), s[1-x]);
			} else {
				beg = p[0]->beg > p[1]->beg? p[0]->beg : p[1]->beg;
				end = tract + beg;
				assert(end != ab->m);
				if (beg == p[0]->beg) {
					if (beg == p[1]->beg) x = (drand48() < 0.5);
					else x = 0;
				} else x = 1;
				coalesce(ab, recombine(ab, end, s[x], 0), s[1-x]);
			}
		}
	}
	{ // print the sequence at the root
		arg_S1_t *p;
		int i;
		k = kh_get(S, ab->S, ab->sid-1); p = &kh_val(ab->S, k);
		assert(p->beg == 0 && p->end == ab->m);
		printf("S\t%d\t", p->nid);
		for (i = p->beg; i < p->end; ++i)
			putchar('0' + p->s[i - p->beg]);
		putchar('\n');
	}
	fprintf(stderr, "[arg_build] # recombinations: %d\n", g_n_recombs);
}

int main_build(int argc, char *argv[])
{
	arg_data_t *ad;
	arg_build_t *ab;
	int c, flag = 0;
	long seed = 11;
	double ratio = 0.0;
	while ((c = getopt(argc, argv, "qds:br:L")) >= 0) {
		switch (c) {
		case 'd': flag |= AB_FSHOW_SID; break;
		case 'q': flag |= AB_FSHOW_SEQ; break;
		case 'b': flag |= AB_FBACKWARD; break;
		case 'L': flag |= AB_FLONGEST;  break;
		case 'r': ratio = atof(optarg); break;
		case 's':
			seed = atol(optarg);
			if (seed <= 0) seed = getpid() ^ time(0);
			break;
		default:
			fprintf(stderr, "[main_arg_build] unrecognized option.\n");
			return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fastARG build [options] <seq.txt>\n\n");
		fprintf(stderr, "Options: -s INT    seed [%ld]\n", seed);
		fprintf(stderr, "         -r FLOAT  fraction of tracts [0.0]\n");
		fprintf(stderr, "         -L        prioritize on length rather than sum of lengths\n");
		fprintf(stderr, "         -d        print sid (for debugging)\n");
		fprintf(stderr, "         -q        print sequence (for debugging)\n");
		fprintf(stderr, "         -b        right-to-left construction\n");
		fprintf(stderr, "\n");
		return 1;
	}
	srand48(seed);
	printf("E\t%ld\n", seed);
	ad = arg_data_read(argv[optind]);
	ab = init_build(ad, flag, ratio);
	arg_data_destroy(ad);
	arg_build(ab);
	destroy_build(ab);
	return 0;
}
