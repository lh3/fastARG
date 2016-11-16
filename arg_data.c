#include <zlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "kseq.h"
#include "arg_data.h"

KSTREAM_INIT(gzFile, gzread, 4096)

arg_data_t *arg_data_read(const char *fn)
{
	kstream_t *ks;
	kstring_t *str;
	gzFile fp;
	arg_data_t *ad;
	int max_m, dret, i;

	ad = (arg_data_t*)calloc(1, sizeof(arg_data_t));
	ad->is_phased = 1;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);

	max_m = 0;
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int pos = atoi(str->s);
		int n[4];
		ks_getuntil(ks, 0, str, &dret);
		if (ad->n == 0) {
			ad->n = strlen(str->s);
			ad->seq = (char**)calloc(ad->n, sizeof(char*));
		} else if (ad->n != strlen(str->s)) {
			fprintf(stderr, "[arg_data_read] wrong input: variable individuals?\n");
			exit(1);
		}
		// change scale
		n[0] = n[1] = n[2] = n[3] = 0;
		for (i = 0; i < ad->n; ++i) {
			str->s[i] = (str->s[i] < '0' || str->s[i] > '2')? 3 : str->s[i] - '0';
			++n[(int)str->s[i]];
		}
		if (n[2] > 0) ad->is_phased = 0;
		if (ad->is_phased && n[2]%2 != 0) {
			fprintf(stderr, "[arg_data_read] wrong input: unphased input with odd number of haplotypes\n");
			exit(1);
		}
		if (n[3] == ad->n)
			for (i = 0; i < ad->n; ++i) str->s[i] = 0;
		if (n[2] > 0) {
			for (i = 0; i < ad->n; i += 2) { // check
				if ((str->s[i] == 2 && str->s[i+1] != 2) || (str->s[i] != 2 && str->s[i+1] == 2)) {
					if (str->s[i] == 2) str->s[i+1] = 2;
					else str->s[i] = 2;
					fprintf(stderr, "[arg_read_data] fix inconsistent unphased character at (%d,%d)\n", pos, i);
				}
			}
			if (n[2] + n[3] == ad->n) { // IMPORTANT: arbitrary phasing the first unphased character
				for (i = 0; i < ad->n; i += 2) {
					if (str->s[i] == 2) {
						str->s[i] = 0; str->s[i+1] = 1;
						break;
					}
				}
			}
		}

		if (ad->m == max_m) {
			max_m = max_m? max_m << 1 : 16;
			ad->pos = (int*)realloc(ad->pos, sizeof(int) * max_m);
			for (i = 0; i < ad->n; ++i)
				ad->seq[i] = (char*)realloc(ad->seq[i], max_m);
		}
		ad->pos[ad->m] = pos;
		for (i = 0; i < ad->n; ++i) ad->seq[i][ad->m] = str->s[i]; // copy over
		++ad->m;
	}

	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	fprintf(stderr, "[arg_data_read] %d samples; %d sites\n", ad->n, ad->m);
	if (!ad->is_phased) fprintf(stderr, "[arg_data_read] the input is unphased\n");
	return ad;
}

void arg_data_print(const arg_data_t *ad)
{
	int i, j;
	for (i = 0; i < ad->m; ++i) {
		printf("%d\t", ad->pos[i]);
		for (j = 0; j < ad->n; ++j)
			putchar(ad->seq[j][i] + '0');
		putchar('\n');
	}
}

void arg_data_destroy(arg_data_t *ad)
{
	int i;
	if (ad == 0) return;
	for (i = 0; i < ad->n; ++i) free(ad->seq[i]);
	free(ad->seq); free(ad->pos);
	free(ad);
}

/* Requirement:
	1. seq[2*k] and seq[2*k+1] constitute a diploid sequence
	2. if seq[0]!=seq[1], seq[2*k]!=seq[2*k+1] always stands
*/

char *arg_minimize_switches(int m, int n, char **seq)
{
	int i, j, k, n_het = 0, *f, *bt;
	char *last_char, *s, *ret;
	assert(n%2 == 0);
	// count # hets
	for (i = k = 0; i < m; ++i)
		if (seq[0][i] != seq[1][i]) ++k;
	n_het = k;
	// DP
	f = (int*)calloc(n_het * 2, sizeof(int));
	bt = (int*)calloc(n_het * 2, sizeof(int));
	last_char = (char*)calloc(n/2, 1);
	for (i = k = 0; i < m; ++i) {
		int cnt[2];
		if (seq[0][i] == seq[1][i]) continue; // not a het
		cnt[0] = cnt[1] = 0;
		for (j = 0; j < n; j += 2) {
			++cnt[(last_char[j/2] == seq[j][i])? 0 : 1];
			last_char[j/2] = seq[j][i];
		}
		if (k > 0) { // not the first het
			int y[2];
			y[0] = f[(k-1)<<1|0] + cnt[0];
			y[1] = f[(k-1)<<1|1] + cnt[1];
			f[k<<1|0] = y[0] > y[1]? y[0] : y[1];
			bt[k<<1|0] = y[0] > y[1]? 0 : 1;
			y[0] = f[(k-1)<<1|0] + cnt[1];
			y[1] = f[(k-1)<<1|1] + cnt[0];
			f[k<<1|1] = y[0] > y[1]? y[0] : y[1];
			bt[k<<1|1] = y[0] > y[1]? 0 : 1;
		} else f[0] = 0, f[1] = -0x7fffffff;
		++k;
	}
	free(last_char);
	// backtrace
	s = (char*)calloc(n_het, 1);
	j = f[(n_het-1)<<1|0] > f[(n_het-1)<<1|1]? 0 : 1;
	for (k = n_het - 1; k >= 0; --k) {
		s[k] = j;
		j = bt[k<<1|j];
	}
	free(bt); free(f);
	// write ret seq
	ret = (char*)calloc(m, 1);
	for (i = k = 0; i < m; ++i)
		ret[i] = (seq[0][i] == seq[1][i])? seq[0][i] : s[k++];
	free(s);
	return ret;
}

int main_fadmerge(int argc, char *argv[])
{
	arg_data_t *ret, **ad;
	int n_fad, i, j;
	char **seq, *s;
	if (argc < 4) {
		fprintf(stderr, "Usage: fastARG fadmerge <in1.fad> <in2.fad> <in3.fad> ...\n");
		return 1;
	}
	// load fads
	n_fad = argc - 1;
	ad = (arg_data_t**)calloc(n_fad, sizeof(void*));
	for (i = 0; i < n_fad; ++i)
		ad[i] = arg_data_read(argv[i + 1]);
	// merge
	ret = (arg_data_t*)calloc(1, sizeof(arg_data_t));
	ret->n = ad[0]->n; ret->m = ad[0]->m;
	ret->pos = (int*)calloc(ret->m, sizeof(int));
	ret->seq = (char**)calloc(ret->n, sizeof(void*));
	seq = (char**)calloc(n_fad * 2, sizeof(void*));
	for (j = 0; j < ret->n; j += 2) {
		for (i = 0; i < n_fad; ++i) {
			seq[2*i] = ad[i]->seq[j];
			seq[2*i+1] = ad[i]->seq[j+1];
		}
		s = arg_minimize_switches(ret->m, n_fad * 2, seq);
		ret->seq[j] = s;
		s = ret->seq[j+1] = (char*)calloc(ret->m, 1);
		for (i = 0; i < ret->m; ++i)
			s[i] = (seq[0][i] == seq[1][i])? seq[0][i] : 1 - ret->seq[j][i];
	}
	arg_data_print(ret);
	// destroy
	for (i = 0; i < n_fad; ++i) arg_data_destroy(ad[i]);
	free(ad); free(seq);
	arg_data_destroy(ret);
	return 0;
}

// sample m integers from [0..n-1] without replacement; ret is of size m
void sample_no_rep(int n, int m, int *ret)
{
	int i, j, k;
	assert(m <= n);
	for (i = n, j = m, k = 0; j > 0; --j) {
		double p = 1.0, x = drand48();
		while (x < p) p -= p * j / (i--);
		ret[k++] = n - i - 1;
	}
}

int main_mixphase(int argc, char *argv[])
{
	int c, n_fuzzy = 0, *inds = 0, i, j;
	double r_fuzzy = 0.5;
	arg_data_t *ad;
	srand48(time(0) ^ getpid());
	while ((c = getopt(argc, argv, "n:")) >= 0) {
		switch (c) {
		case 'n':
			r_fuzzy = atof(optarg);
			if (r_fuzzy >= 1.0) n_fuzzy = (int)(r_fuzzy + .499);
			break;
		default: return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fastARG mixphase [-n ratio|n] <in.fad>\n");
		return 1;
	}
	ad = arg_data_read(argv[optind]);
	if (n_fuzzy > ad->n/2) {
		fprintf(stderr, "[main_mixphase] only %d individuals are in the input.\n", ad->n/2);
		return 1;
	}
	if (n_fuzzy == 0) n_fuzzy = (int)(ad->n/2 * r_fuzzy + .499);
	if (inds == 0) {
		inds = (int*)calloc(n_fuzzy, sizeof(int));
		sample_no_rep(ad->n/2, n_fuzzy, inds);
	}
	// mixphase: core loop
	for (i = 0; i < n_fuzzy; ++i) {
		char *s[2];
		s[0] = ad->seq[inds[i]*2];
		s[1] = ad->seq[inds[i]*2 + 1];
		for (j = 0; j < ad->m; ++j)
			if (s[0][j] != s[1][j] && s[0][j] < 2 && s[1][j] < 2)
				s[0][j] = s[1][j] = 2;
	}
	arg_data_print(ad);
	// free
	free(inds);
	arg_data_destroy(ad);
	return 0;
}
