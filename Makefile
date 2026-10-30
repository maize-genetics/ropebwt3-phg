CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O3
CPPFLAGS=
INCLUDES=
OBJS=		libsais.o libsais64.o kalloc.o kthread.o misc.o io.o rld0.o bre.o rle.o rope.o mrope.o \
			dawg.o fm-index.o ssa.o sais-ss.o build.o search.o bwa-sw.o lift.o ps4g.o hitcount.o pairing.o isa.o extract.o snap.o errprof.o edstats.o
PROG=		ropebwt3
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(omp),0)
	CPPFLAGS=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
	LIBS+=-fopenmp
endif

.SUFFIXES:.c .o
.PHONY:all clean depend test

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ropebwt3:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

rld0.o:rld0.c rld0.h bre.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DRLD_HAVE_BRE $(INCLUDES) $< -o $@

test/test_ps4g.o:test/test_ps4g.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

# ps4g.o now calls rb3_lift_nearest_ref/rb3_lift_ternary_state (lift.o), even though
# this test only ever exercises them with ridx=0 -- the linker still needs the
# symbols. lift.o pulls in the FM-index/SSA/kalloc/kthread objects transitively
# (see test/test_lift's own link line below).
test/test_ps4g:test/test_ps4g.o ps4g.o misc.o lift.o fm-index.o ssa.o kalloc.o kthread.o io.o rld0.o bre.o rle.o rope.o mrope.o
		$(CC) $(CFLAGS) $^ -o $@ -lpthread -lz -lm

test/test_hitcount.o:test/test_hitcount.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_hitcount:test/test_hitcount.o hitcount.o
		$(CC) $(CFLAGS) $^ -o $@ -lpthread

test/test_lift.o:test/test_lift.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_lift:test/test_lift.o lift.o ps4g.o fm-index.o ssa.o kalloc.o kthread.o misc.o io.o rld0.o bre.o rle.o rope.o mrope.o
		$(CC) $(CFLAGS) $^ -o $@ -lpthread -lz -lm

test/test_pairing.o:test/test_pairing.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_pairing:test/test_pairing.o pairing.o
		$(CC) $(CFLAGS) $^ -o $@ -lm

test/test_extract.o:test/test_extract.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_extract:test/test_extract.o snap.o isa.o lift.o ps4g.o fm-index.o ssa.o kalloc.o kthread.o misc.o io.o rld0.o bre.o rle.o rope.o mrope.o sais-ss.o libsais.o libsais64.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

test/test_errprof.o:test/test_errprof.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_errprof:test/test_errprof.o errprof.o
		$(CC) $(CFLAGS) $^ -o $@ -lm

test/test_edstats.o:test/test_edstats.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -I. $< -o $@

test/test_edstats:test/test_edstats.o edstats.o errprof.o
		$(CC) $(CFLAGS) $^ -o $@ -lm

test:$(PROG) test/test_ps4g test/test_hitcount test/test_lift test/test_pairing test/test_extract test/test_errprof test/test_edstats
		@mkdir -p test/output
		./test/test_ps4g test/output
		./test/test_hitcount
		./test/test_lift
		./test/test_pairing
		./test/test_extract test/output
		./test/test_errprof
		./test/test_edstats
		cd test && ./run_integration.sh

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM test/*.o test/test_ps4g test/test_hitcount test/test_lift test/test_pairing test/test_extract test/test_errprof test/test_edstats test/output test/tmp

# Regenerate the header dependencies below (everything after "DO NOT DELETE") with the
# compiler's -MM; run `make depend` after adding a source file or changing an #include.
depend:
		(LC_ALL=C; export LC_ALL; sed '/^# DO NOT DELETE/q' Makefile > Makefile.dep && echo >> Makefile.dep && \
		 for f in *.c test/*.c; do $(CC) -MM -MT $${f%.c}.o -I. $(CPPFLAGS) -DRLD_HAVE_BRE $$f || exit 1; done >> Makefile.dep && \
		 mv Makefile.dep Makefile)

# DO NOT DELETE

bre.o: bre.c bre.h
build.o: build.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h rle.h \
 bre.h ketopt.h kthread.h
bwa-sw.o: bwa-sw.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h \
 align.h kalloc.h dawg.h khashl-km.h ksort.h
dawg.o: dawg.c dawg.h kalloc.h libsais.h io.h rb3priv.h khashl-km.h
edstats.o: edstats.c edstats.h errprof.h
errprof.o: errprof.c errprof.h
extract.o: extract.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h \
 ketopt.h kalloc.h kthread.h lift.h ps4g.h isa.h snap.h
fm-index.o: fm-index.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h \
 rle.h kthread.h kalloc.h khashl-km.h
hitcount.o: hitcount.c hitcount.h
io.o: io.c rb3priv.h io.h kseq.h
isa.o: isa.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h ketopt.h \
 isa.h
kalloc.o: kalloc.c kalloc.h
kthread.o: kthread.c kthread.h
libsais.o: libsais.c libsais.h
libsais64.o: libsais64.c libsais.h libsais64.h
lift.o: lift.c fm-index.h rb3priv.h rld0.h mrope.h rope.h io.h ketopt.h \
 kalloc.h kthread.h lift.h
main.o: main.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h ketopt.h
misc.o: misc.c rb3priv.h
mrope.o: mrope.c mrope.h rope.h rle.h
pairing.o: pairing.c pairing.h
ps4g.o: ps4g.c rb3priv.h ps4g.h io.h lift.h kseq.h
rld0.o: rld0.c rld0.h bre.h
rle.o: rle.c rle.h
rope.o: rope.c rle.h rope.h
sais-ss.o: sais-ss.c rb3priv.h libsais.h libsais64.h
search.o: search.c fm-index.h rb3priv.h rld0.h mrope.h rope.h io.h \
 align.h ketopt.h kthread.h kalloc.h lift.h ps4g.h hitcount.h pairing.h \
 errprof.h edstats.h khashl-km.h
snap.o: snap.c rb3priv.h snap.h
ssa.o: ssa.c rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h kalloc.h \
 kthread.h ketopt.h ksort.h
test/test_edstats.o: test/test_edstats.c edstats.h
test/test_errprof.o: test/test_errprof.c errprof.h
test/test_extract.o: test/test_extract.c rb3priv.h fm-index.h rb3priv.h \
 rld0.h mrope.h rope.h io.h io.h kalloc.h ps4g.h lift.h lift.h isa.h \
 fm-index.h snap.h
test/test_hitcount.o: test/test_hitcount.c hitcount.h
test/test_lift.o: test/test_lift.c rb3priv.h io.h rb3priv.h ps4g.h io.h \
 lift.h lift.h
test/test_pairing.o: test/test_pairing.c pairing.h
test/test_ps4g.o: test/test_ps4g.c rb3priv.h ps4g.h io.h rb3priv.h lift.h
