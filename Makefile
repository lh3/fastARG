CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -m64
CXXFLAGS=	$(CFLAGS)
DFLAGS=		#-DNDEBUG
OBJS=		arg_data.o arg_build.o arg.o arg_check.o
PROG=		fastARG
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

fastARG:$(OBJS) arg_main.o
		$(CC) $(CFLAGS) $(OBJS) arg_main.o -o $@ -lm -lz

arg_data.o:kseq.h arg_data.h
arg.o:kseq.h
arg_build.o:khash.h ksort.h kvec.h arg_data.h
arg.o:arg.h arg_data.h

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
