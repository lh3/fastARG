#ifndef LH3_ARG_DATA_H
#define LH3_ARG_DATA_H

typedef struct {
	int n, m; // n haplotypes, m loci
	int *pos;
	int is_phased;
	char **seq; // IMPORTANT: 0->hom1, 1->hom2, 2->het, 3->N
} arg_data_t;

arg_data_t *arg_data_read(const char *fn);
void arg_data_print(const arg_data_t *ad);
void arg_data_destroy(arg_data_t *ad);

#endif
