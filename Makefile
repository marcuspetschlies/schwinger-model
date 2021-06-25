# ##########
work=.
srcdir=..
gitdir=/home/marcuspe/notes/vl-lqcd-ss2021/exercises

CXX=g++ -fopenmp
CXXFLAGS= -Wall -pedantic -g -O3 -DHAVE_OPENMP -DTIMER=1

CCDEP = g++
DEPFLAGS = -MM

INCLUDE = -I$(srcdir) 

LIBS = -lm

LDFLAGS = -lm 

LINK = $(CXX) -o $@ ${LDFLAGS}
COMPILE = ${CXX} $(INCLUDE) -o $@ ${CXXFLAGS}

MODULES = ranlxd operators solver set_default gauge meas

HEADERS = ranlxd table_init_d table_init_i table_init_u operators solver set_default gauge meas

PROGRAM = schwinger

all: dep $(PROGRAM) 


# ##########

$(addsuffix .d,$(MODULES)): %.d: ${srcdir}/%.c
	 @ $(CCDEP) ${DEPFLAGS} ${INCLUDE} $< > $@

$(addsuffix .d,$(PROGRAM)): %.d: ${srcdir}/%.c
	 @ $(CCDEP) ${DEPFLAGS} ${INCLUDE} $< > $@

dep: $(addsuffix .d,$(MODULES) ${PROGRAM})

$(addsuffix .o,${MODULES}): %.o: ${srcdir}/%.c $(addprefix ${srcdir}/, $(addsuffix .h, ${HEADERS})) %.d
	${COMPILE} ${OPTARGS} -c $< 

$(addsuffix .o,${PROGRAM}): %.o: ${srcdir}/%.c %.d
	${COMPILE} ${OPTARGS} -c $< 

${PROGRAM}: %: %.o gitversion.c $(addsuffix .o,${MODULES})
	${LINK}  $(addsuffix .o,${MODULES}) $@.o $(LIBS)

gitversion.c: ${gitdir}/.git/HEAD ${gitdir}/.git/index
	echo "namespace cvc { const char *g_gitversion = \"$(shell cd ${srcdir} && git rev-parse HEAD && cd -)\"; }" > ${srcdir}/gitversion.c

# ##########


clean:
	rm -f *~ *.o *.d $(PROGRAM) ${srcdir}/gitversion.c

.PHONY: clean

# ##########
