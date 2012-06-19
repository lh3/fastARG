#include <stdio.h>
#include <string.h>
#include "arg_main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.4.0-31 (r1264)"
#endif

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: fastARG\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   fastARG <command> <arguments>\n\n");
	fprintf(stderr, "Command: build         construct ARG\n");
	fprintf(stderr, "         leafseq       get haplotypes from ARG\n");
	fprintf(stderr, "         chkphase      check phasing accuracy\n");
	fprintf(stderr, "         fadmerge      merge multiple fad files\n");
	fprintf(stderr, "         regtest       regression test\n");
	fprintf(stderr, "         relabel       relabel nodes\n");
	fprintf(stderr, "         mixphase      mix phase\n");
	fprintf(stderr, "         del           delete nodes\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "build") == 0) return main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "leafseq") == 0) return main_leafseq(argc-1, argv+1);
	else if (strcmp(argv[1], "chkphase") == 0) return main_chkphase(argc-1, argv+1);
	else if (strcmp(argv[1], "fadmerge") == 0) return main_fadmerge(argc-1, argv+1);
	else if (strcmp(argv[1], "mixphase") == 0) return main_mixphase(argc-1, argv+1);
	else if (strcmp(argv[1], "regtest") == 0) return main_regtest(argc-1, argv+1);
	else if (strcmp(argv[1], "del") == 0) return main_del(argc-1, argv+1);
	else if (strcmp(argv[1], "relabel") == 0) return main_relabel(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
