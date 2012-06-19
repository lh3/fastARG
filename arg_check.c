#include <stdio.h>
#include <stdlib.h>
#include "arg_data.h"

static void chkphase(const arg_data_t *src, arg_data_t *dst)
{
	int i, j;
	int n_switch = 0, n_het = 0;
	if (src->m != dst->m || src->n != dst->n) {
		fprintf(stderr, "[chkphase] inconsistent src and dst.\n");
		exit(1);
	}
	if (src->n%2 != 0) {
		fprintf(stderr, "[chkphase] # haplotypes is not even.\n");
		exit(1);
	}
	for (j = 0; j < src->n; j += 2) {
		char *ss0 = src->seq[j], *ss1 = src->seq[j+1];
		char *sd0 = dst->seq[j], *sd1 = dst->seq[j+1];
		int last = -1, curr, x = 0;
		for (i = 0; i < src->m; ++i) {
			if (ss0[i] == ss1[i]) continue;
			if (sd0[i] == sd1[i]) {
				fprintf(stderr, "[chkphase] inconsistent phase at (%d,%d)\n", j, i);
				exit(1);
			}
			++n_het; ++x;
			curr = (sd0[i] == ss0[i])? 0 : 1;
			if (last >= 0 && last != curr) {
				++n_switch;
				x = 1;
			}
			last = curr;
		}
	}
	fprintf(stderr, "[chkphase] # heterozygotes: %d\n", n_het);
	fprintf(stderr, "[chkphase] # switches:      %d\n", n_switch);
	fprintf(stderr, "[chkphase] switch error:    %.3lf\n", (double)n_switch/(n_het - src->n/2));
}

int main_chkphase(int argc, char *argv[])
{
	arg_data_t *src, *dst;
	if (argc < 3) {
		fprintf(stderr, "Usage: fastARG chkphase <src.dat> <dst.dat>\n");
		return 1;
	}
	src = arg_data_read(argv[1]);
	dst = arg_data_read(argv[2]);
	chkphase(src, dst);
	arg_data_destroy(src);
	arg_data_destroy(dst);
	return 0;
}
