#ifndef LH3_ARG_H
#define LH3_ARG_H

#include <stdint.h>

typedef struct {
	/* if n_nei==0, the node is disconnected. */
	int nid, n_nei;
	/* for a C node, nei[0,1] are the children and nei[2] is the parent;
	 * for an R node, nei[0] is the child and nei[1,2] are the parents;
	 * for a leaf node, nei[0] is the parent. */
	int nei[3];
	/* x==0 for a C node and beg<x<end for an R node */
	int beg, x, end;
	/* mutations on incoming arcs */
	int n_mut, *mut;
} argnode_t;

typedef struct {
	/* m=#sites; n=#haplotypes */
	int max_size, n, m, root;
	uint64_t *rootseq;
	argnode_t *node;
} arg_t;

#define arg_getseq(s, i) (s[(i)>>6]>>(~(i)&0x3f)&1)
#define arg_setseq0(s, i) (s[(i)>>6] &= ~(1ull<<(~(i)&0x3f)))
#define arg_setseq1(s, i) (s[(i)>>6] |= 1ull<<(~(i)&0x3f))
#define arg_xorseq(s, i) (s[(i)>>6] ^= 1ull<<(~(i)&0x3f))
#define arg_del1(p) do { free((p)->mut); (p)->mut = 0; (p)->n_mut = (p)->n_nei = 0; } while (0)
#define arg_childreg(p, q, _b, _e) do {						\
		*(_b) = (q)->beg; *(_e) = (q)->end;					\
		if ((q)->x) {										\
			if ((q)->nei[1] == (p)->nid) *(_e) = (q)->x;	\
			else *(_b) = (q)->x;							\
		}													\
	} while (0)

#ifdef __cplusplus
extern "C" {
#endif

	arg_t *arg_init();
	arg_t *arg_dup(const arg_t *a);
	void arg_destroy(arg_t *a);

	arg_t *arg_load(gzFile fp);
	void arg_print(const arg_t *a);
	int arg_check(const arg_t *a);

	uint64_t **arg_leaf_seq(const arg_t *a);
	int arg_del(arg_t *a, int leaf_id);

	void arg_sort(arg_t *a);

#ifdef __cplusplus
}
#endif

#endif
