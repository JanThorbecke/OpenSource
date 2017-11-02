# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
LIBS    += -L$L -lgenfft -lm $(LIBSM)
#LIBS    += -L$L -lgenfft -lm -lc
OPTC = -g -Wall -fsignaling-nans -O0
#OPTC += -fopenmp -Waddress
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#LDFLAGS += -Mprof=lines

all: raytime 

PRG = raytime

SRCC	= $(PRG).c \
		JespersRayTracer.c \
		getParameters.c  \
		getWaveletInfo.c  \
		writeSrcRecPos.c  \
		readModel.c  \
		getWaveletHeaders.c  \
		verbosepkg.c  \
        getModelInfo.c  \
		wallclock_time.c  \
        recvPar.c  \
        writesufile.c  \
        gaussGen.c  \
        name_ext.c  \
		atopkge.c \
		docpkge.c \
		threadAffinity.c \
		getpars.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) raytime.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o raytime $(OBJC) $(LIBS)

install: raytime 
	cp raytime $B

clean:
		rm -f core $(OBJC) $(OBJM) raytime 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/raytime 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



