# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
LIBS    += -L$L -lgenfft $(LIBSM)
#LIBS    += -L$L -lgenfft -lm -lc
#OPTC = -g -Wall -fsignaling-nans -O0 -fopenmp
#OPTC += -g3 
#OPTC += -fopenmp -Waddress
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#OPTC += -qopt-report
#LDFLAGS += -Mprof=lines

all: fdelmodc3D

PRG = fdelmodc3D

SRCC	= $(PRG).c \
		acoustic2_3D.c \
		acoustic4_3D.c \
		acoustic6_3D.c \
		acousticSH4_3D.c \
		acoustic4_qr_3D.c \
		elastic4dc_3D.c \
		elastic4_3D.c \
		defineSource3D.c  \
		getParameters3D.c  \
		getWaveletInfo3D.c  \
		getModelInfo3D.c  \
		applySource3D.c  \
		getRecTimes3D.c  \
		getBeamTimes3D.c  \
		writeSnapTimes3D.c  \
		writeRec3D.c  \
		writeSrcRecPos3D.c  \
		decomposition.c  \
		fileOpen.c  \
		recvPar3D.c  \
		readModel3D.c  \
		sourceOnSurface3D.c  \
		getWaveletHeaders3D.c  \
		boundaries3D.c  \
		verbosepkg.c  \
		writesufile3D.c  \
		writesufile.c  \
		gaussGen.c  \
		spline3.c  \
		CMWC4096.c  \
		wallclock_time.c  \
		name_ext.c  \
		atopkge.c \
		docpkge.c \
		threadAffinity.c \
		alloc3D.c \
		getpars.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) fdelmodc3D.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o fdelmodc3D $(OBJC) $(LIBS)

install: fdelmodc3D 
	cp fdelmodc3D $B

clean:
		rm -f core $(OBJC) $(OBJM) fdelmodc3D

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/fdelmodc3D


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



