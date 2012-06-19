#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include "arg.h"
#include "kvec.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 4096)

arg_t *arg_init()
{
	arg_t *a = (arg_t*)calloc(1, sizeof(arg_t));
	return a;
}

void arg_destroy(arg_t *a)
{
	int i;
	for (i = 0; i <= a->root; ++i)
		free(a->node[i].mut);
	free(a->node); free(a->rootseq);
	free(a);
}

arg_t *arg_load(gzFile fp)
{
	kstream_t *ks;
	kstring_t *str;
	int dret, lineno = 0;
	arg_t *a;
	a = arg_init();
	ks = ks_init(fp);
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int n1, n2, i;
		++lineno;
		if (str->l != 1) {
			fprintf(stderr, "[arg_load] invalid initial character at line %d\n", lineno);
			exit(1);
		}
		if (str->s[0] == 'C' || str->s[0] == 'R') {
			argnode_t *an1, *an2;
			int beg, end, n_mut;
			// read
			ks_getuntil(ks, 0, str, &dret); n1 = atoi(str->s);
			ks_getuntil(ks, 0, str, &dret); n2 = atoi(str->s);
			if (n1 <= n2) {
				fprintf(stderr, "[arg_load] invalid edge (%d)\n", lineno);
				exit(1);
			}
			// add to the node array
			if (n1 + 1 > a->max_size) {
				int old_size = a->max_size;
				a->max_size = n1 + 1;
				kroundup32(a->max_size);
				a->node = (argnode_t*)realloc(a->node, sizeof(argnode_t) * a->max_size);
				memset(a->node + old_size, 0, (a->max_size - old_size) * sizeof(argnode_t));
			}
			an1 = a->node + n1; an2 = a->node + n2;
			an1->nid = n1; an2->nid = n2;
			if (an1->n_nei == 3 || an2->n_nei == 3) {
				fprintf(stderr, "[arg_load] multifurcated node: %d or %d (%d)\n", n1, n2, lineno);
				exit(1);
			}
			an1->nei[an1->n_nei++] = n2; an2->nei[an2->n_nei++] = n1;
			// read intervals
			ks_getuntil(ks, 0, str, &dret); beg = atoi(str->s);
			ks_getuntil(ks, 0, str, &dret); end = atoi(str->s);
			if (an2->end) { // a recombination node
				if (an2->end == beg) an2->x = an2->end, an2->end = end;
				else if (an2->beg == end) { // and also swap
					int x = an2->nei[1]; an2->nei[1] = an2->nei[2]; an2->nei[2] = x;
					an2->x = an2->beg, an2->beg = beg;
				} else {
					fprintf(stderr, "[arg_load] inconsisten interval at node %d (%d)\n", n2, lineno);
					exit(1);
				}
			} else an2->beg = beg, an2->end = end;
			// read mutations
			ks_getuntil(ks, 0, str, &dret); n_mut = atoi(str->s);
			if (n_mut) {
				an2->mut = (int*)realloc(an2->mut, sizeof(int) * (an2->n_mut + n_mut));
				for (i = 0; i < n_mut; ++i) {
					ks_getuntil(ks, 0, str, &dret);
					an2->mut[an2->n_mut++] = atoi(str->s);
				}
			}
			// save interval and swap if necessary
		} else if (str->s[0] == 'N') {
			ks_getuntil(ks, 0, str, &dret); a->n = atoi(str->s);
			ks_getuntil(ks, 0, str, &dret); a->m = atoi(str->s);
		} else if (str->s[0] == 'S') {
			ks_getuntil(ks, 0, str, &dret); a->root = atoi(str->s);
			ks_getuntil(ks, 0, str, &dret);
			if (str->l != a->m) {
				fprintf(stderr, "[arg_load] inconsistent root sequence (%d)\n", lineno);
				exit(1);
			}
			a->rootseq = (uint64_t*)calloc((a->m + 63) / 64, 8);
			for (i = 0; i < a->m; ++i)
				if (str->s[i] == '1') arg_setseq1(a->rootseq, i);
		} else ks_getuntil(ks, '\n', str, &dret);
	}
	free(str->s); free(str);
	ks_destroy(ks);
	return a;
}

/* q is the child, p is the parent, s is the sequence at q, t is the
 * sequence at p, and nid is the node id of p */
static inline uint64_t *propagate_seq(int nid, argnode_t *p, argnode_t *q, uint64_t *s, uint64_t *t)
{
	int i, beg, end;
	beg = q->beg; end = q->end;
	if (q->x) { // then q cannot be a leaf
		if (nid == q->nei[1]) end = q->x;
		else beg = q->x;
	}
	if (t == 0) t = (uint64_t*)calloc((q->end - q->beg + 63) / 64, 8);
	for (i = beg; i < end; ++i)
		if (arg_getseq(s, i - p->beg))
			arg_setseq1(t, i - q->beg);
	for (i = 0; i < q->n_mut; ++i)
		if (q->mut[i] >= beg && q->mut[i] < end)
			arg_xorseq(t, q->mut[i] - q->beg);
	return t;
}

uint64_t **arg_leaf_seq(const arg_t *a)
{
	uint64_t **seq;
	int j;
	seq = (uint64_t**)calloc(a->root + 1, sizeof(void*));
	seq[a->root] = (uint64_t*)calloc((a->m + 63) / 64, 8);
	memcpy(seq[a->root], a->rootseq, (a->m + 63) / 64 * 8);
	for (j = a->root; j >= a->n; --j) {
		argnode_t *p = a->node + j;
		if (p->n_nei == 0) continue;
		seq[p->nei[0]] = propagate_seq(j, p, a->node + p->nei[0], seq[j], seq[p->nei[0]]);
		if (p->x == 0) seq[p->nei[1]] = propagate_seq(j, p, a->node + p->nei[1], seq[j], seq[p->nei[1]]);
		free(seq[j]);
		seq[j] = 0;
	}
	seq = (uint64_t**)realloc(seq, sizeof(void*) * a->n);
	return seq;
}

int main_leafseq(int argc, char *argv[])
{
	arg_t *a;
	gzFile fp;
	uint64_t **seq;
	int i, j;
	if (argc == 1) {
		fprintf(stderr, "Usage: fastARG leafseq <in.arg>\n");
		return 1;
	}
	fp = strcmp(argv[1], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	a = arg_load(fp);
	gzclose(fp);
	seq = arg_leaf_seq(a);
	for (i = 0; i < a->m; ++i) {
		printf("%d\t", i);
		for (j = 0; j < a->n; ++j)
			if (seq[j]) putchar(arg_getseq(seq[j], i) + '0');
		putchar('\n');
	}
	// FIXME: memory leak
	arg_destroy(a);
	return 0;
}

int arg_check(const arg_t *a)
{
	int i;
	for (i = 0; i < a->root; ++i) {
		argnode_t *q, *p = a->node + i;
		if (p->n_nei == 0) continue;
		if (p->x == 0) {
			q = a->node + p->nei[p->n_nei-1];
			if (i != q->nei[0] && i != q->nei[1]) break;
		} else {
			q = a->node + p->nei[1];
			if (i != q->nei[0] && i != q->nei[1]) break;
			q = a->node + p->nei[2];
			if (i != q->nei[0] && i != q->nei[1]) break;
		}
	}
	if (i == a->root) return 0;
	fprintf(stderr, "%d\n", i);
	return -1;
}

void arg_print(const arg_t *a)
{
	int i, j;
	printf("N\t%d\t%d\n", a->n, a->m);
	for (i = 0; i < a->root; ++i) {
		argnode_t *p = a->node + i;
		if (p->n_nei == 0) continue;
		if (p->x == 0) { // one parent
			putchar(a->node[p->nei[p->n_nei-1]].x? 'R' : 'C');
			printf("\t%d\t%d\t%d\t%d\t%d", p->nei[p->n_nei-1], i, p->beg, p->end, p->n_mut);
			for (j = 0; j < p->n_mut; ++j) printf("\t%d", p->mut[j]);
			putchar('\n');
		} else { // two parents
			int n_mut;
			putchar(a->node[p->nei[1]].x? 'R' : 'C');
			printf("\t%d\t%d\t%d\t%d", p->nei[1], i, p->beg, p->x);
			for (j = n_mut = 0; j < p->n_mut; ++j)
				if (p->mut[j] >= p->beg && p->mut[j] < p->x) ++n_mut;
			printf("\t%d", n_mut);
			if (n_mut) {
				for (j = 0; j < p->n_mut; ++j)
					if (p->mut[j] >= p->beg && p->mut[j] < p->x) printf("\t%d", p->mut[j]);
			}
			putchar('\n');
			putchar(a->node[p->nei[2]].x? 'R' : 'C');
			printf("\t%d\t%d\t%d\t%d", p->nei[2], i, p->x, p->end);
			for (j = n_mut = 0; j < p->n_mut; ++j)
				if (p->mut[j] >= p->x && p->mut[j] < p->end) ++n_mut;
			printf("\t%d", n_mut);
			if (n_mut) {
				for (j = 0; j < p->n_mut; ++j)
					if (p->mut[j] >= p->x && p->mut[j] < p->end) printf("\t%d", p->mut[j]);
			}
			putchar('\n');
		}
	}
	printf("S\t%d\t", a->root);
	for (i = 0; i < a->m; ++i)
		putchar(arg_getseq(a->rootseq, i) + '0');
	putchar('\n');
}

int main_regtest(int argc, char *argv[])
{
	arg_t *a;
	gzFile fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: fastARG regtest <in.arg>\n");
		return 1;
	}
	fp = strcmp(argv[1], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	a = arg_load(fp);
	gzclose(fp);
	arg_print(a);
	arg_destroy(a);
	return 0;
}

arg_t *arg_dup(const arg_t *a)
{
	int i;
	arg_t *aa;
	aa = arg_init();
	*aa = *a;
	aa->max_size = a->root + 1;
	aa->rootseq = calloc((a->m + 63)/64, 8);
	memcpy(aa->rootseq, a->rootseq, (a->m+63)/64*8);
	aa->node = calloc(a->root+1, sizeof(argnode_t));
	memcpy(aa->node, a->node, (a->root+1) * sizeof(argnode_t));
	for (i = 0; i <= a->root; ++i) {
		argnode_t *p = aa->node + i;
		if (p->n_mut > 0) {
			p->mut = malloc(p->n_mut * sizeof(int));
			memcpy(p->mut, a->node[i].mut, p->n_mut * sizeof(int));
		}
	}
	return aa;
}

int arg_del(arg_t *a, int id)
{
	int n = 0;
	argnode_t *p, *q;
	kvec_t(int) stack;
	if (id >= a->n) return -1;
	kv_init(stack);
	q = a->node + a->node[id].nei[0]; // parent
	kv_push(int, stack, a->node[id].nei[0]);
	q->nei[q->nei[0] == id? 0 : 1] = -1;
	arg_del1(a->node + id);
	while (kv_size(stack)) {
		int parent, self = kv_pop(stack), is_del = 1;
		if (self < 0) self = -self, is_del = 0;
		p = a->node + self;
		if (is_del) { // delete p
			if (p->x) { // p has two parents
				q = a->node + p->nei[1];
				q->nei[q->nei[0] == self? 0 : 1] = -1;
				kv_push(int, stack, p->nei[1]);
				q = a->node + p->nei[2];
				q->nei[q->nei[0] == self? 0 : 1] = -1;
				kv_push(int, stack, p->nei[2]);
				arg_del1(p);
				++n;
			} else if (self != a->root) { // p has one parent and is not the root
				int child = p->nei[p->nei[0] >= 0? 0 : 1];
				argnode_t *r = a->node + child;
				kv_push(int, stack, -p->nei[2]);
				parent = p->nei[2];
				q = a->node + parent;
				q->nei[q->nei[0] == self? 0 : 1] = child;
				if (r->n_nei == 1) r->nei[0] = parent;
				else r->nei[r->nei[1] == self? 1 : 2] = parent;
				if (p->n_mut) {
					r->mut = (int*)realloc(r->mut, (r->n_mut + p->n_mut) * sizeof(int));
					memcpy(r->mut + r->n_mut, p->mut, p->n_mut * sizeof(int));
					r->n_mut += p->n_mut;
				}
				arg_del1(p);
			}
		} else { // adjust the active region only
			if (p->x) {
				int beg, end;
				arg_childreg(p, a->node + p->nei[0], &beg, &end);
				if (p->beg > beg || p->end < end) {
					fprintf(stderr, "BUG?!\n");
					continue;
				} else if (p->beg == beg && p->end == end) {
					continue;
				} else if (p->x <= beg || p->x >= end) {
					argnode_t *r;
					if (p->x <= beg) {
						kv_push(int, stack, p->nei[1]);
						p->beg = beg;
						r = a->node + p->nei[1];
						r->nei[r->nei[0] == p->nid? 0 : 1] = -1;
					} else {
						kv_push(int, stack, p->nei[2]);
						p->end = end;
						r = a->node + p->nei[2];
						r->nei[r->nei[0] == p->nid? 0 : 1] = -1;
						p->nei[2] = p->nei[1];
					}
					p->x = 0; p->nei[1] = -1; // fake a C node
					kv_push(int, stack, self);
					continue;
				}
				if (p->beg < beg && beg < p->x) {
					p->beg = beg;
					kv_push(int, stack, -p->nei[1]);
				}
				if (p->end > end && end > p->x) {
					p->end = end;
					kv_push(int, stack, -p->nei[2]);
				}
			} else if (self != a->root) {
				int k, beg, end, b[2], e[2];
				for (k = 0; k < 2; ++k)
					arg_childreg(p, a->node + p->nei[k], &b[k], &e[k]);
				beg = b[0] < b[1]? b[0] : b[1];
				end = e[0] > e[1]? e[0] : e[1];
				if (p->beg != beg || p->end != end) {
					p->beg = beg; p->end = end;
					kv_push(int, stack, -p->nei[2]);
				}
			}
		}
	}
	kv_destroy(stack);
	return n;
}

int main_del(int argc, char *argv[])
{
	int i;
	arg_t *a;
	gzFile fp;
	if (argc < 3) {
		fprintf(stderr, "Usage: fastARG del <in.arg> <id1> [<id2> [...]]\n");
		return 1;
	}
	fp = strcmp(argv[1], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	a = arg_load(fp);
	gzclose(fp);
	for (i = 2; i < argc; ++i)
		arg_del(a, atoi(argv[i]));
	arg_print(a);
	arg_destroy(a);
	return 0;	
}

typedef struct {
	int v, p;
} DFS_stack_t;

/* this is a DFS in the finishing order */
int arg_order(const arg_t *a, int *order)
{
	DFS_stack_t *q, *stack;
	int u, *r;
	char *flag;
	r = order;
	flag = calloc(a->root + 1, 1);
	stack = calloc(a->root + 1, sizeof(DFS_stack_t));
	for (u = 0; u < a->n; ++u)
		if (a->node[u].n_nei) *r++ = u;
	stack->v = a->root;
	stack->p = 0;
	flag[0] = 1; // in the stack
	q = stack;
	while (q >= stack) { // stack is not empty
		while (1) {
			argnode_t *p = a->node + q->v;
			int n_child = p->x? 1 : 2;
			if (q->p >= n_child) break;
			u = p->nei[q->p];
			if (flag[u] == 0) { // not visited
				flag[u] = 1; // in the stack
				if (u >= a->n) {
					++q; // push
					q->v = u;
					q->p = 0;
				}
				continue;
			}
			++q->p;
		}
		*r++ = q->v;
		--q;
	}
	free(stack); free(flag);
	return r - order;
}

void arg_sort(arg_t *a)
{
	int *order, *inv, i, x, old_root = a->root;
	argnode_t *old_node = a->node;
	order = calloc(a->root + 1, sizeof(int));
	inv = calloc(a->root + 1, sizeof(int));
	x = arg_order(a, order);
	a->root = x - 1;
	a->node = calloc(x, sizeof(argnode_t));
	for (i = 0; i < x; ++i)
		inv[order[i]] = i;
	free(order);
	for (i = x = 0; i < a->n; ++i)
		if (old_node[i].n_nei) ++x;
	a->n = x;
	for (i = 0; i <= old_root; ++i) {
		argnode_t *q = old_node + i;
		argnode_t *p;
		int k;
		if (q->n_nei == 0) continue;
		p = a->node + inv[q->nid];
		*p = *q;
		for (k = 0; k < q->n_nei; ++k)
			p->nei[k] = inv[q->nei[k]];
		q->n_mut = 0; q->mut = 0;
	}
	free(old_node);
	free(inv);
}

int main_relabel(int argc, char *argv[])
{
	arg_t *a;
	gzFile fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: fastARG relabel <in.arg>\n");
		return 1;
	}
	fp = strcmp(argv[1], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	a = arg_load(fp);
	gzclose(fp);
	arg_sort(a);
	arg_print(a);
	arg_destroy(a);
	return 0;
}
