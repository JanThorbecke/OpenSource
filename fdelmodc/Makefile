# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
#LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC = -g -Wall -fsignaling-nans -O0 -fopenmp
#OPTC += -fopenmp -Waddress
#OPTC +=  -qopt-report=5 -g
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#OPTC += -qopt-report
#LDFLAGS += -Mprof=lines

all: fdelmodc 

PRG = fdelmodc

SRCC	= $(PRG).c \
		acoustic2.c \
		acoustic4.c \
		acousticSH4.c \
		acoustic4_qr.c \
		acoustic6.c \
		viscoacoustic4.c \
		elastic4.c \
		elastic4dc.c \
		elastic6.c \
		viscoelastic4.c \
		defineSource.c  \
		getParameters.c  \
		getWaveletInfo.c  \
		getModelInfo.c  \
		applySource.c  \
		getRecTimes.c  \
		getBeamTimes.c  \
		writeSnapTimes.c  \
		writeRec.c  \
		writeSrcRecPos.c  \
		decomposition.c  \
		fileOpen.c  \
		recvPar.c  \
		readModel.c  \
		sourceOnSurface.c  \
		getWaveletHeaders.c  \
		boundaries.c  \
		verbosepkg.c  \
		writesufile.c  \
		gaussGen.c  \
		spline3.c  \
		CMWC4096.c  \
		wallclock_time.c  \
		name_ext.c  \
		atopkge.c \
		docpkge.c \
		threadAffinity.c \
		getpars.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) fdelmodc.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o fdelmodc $(OBJC) $(LIBS)

install: fdelmodc 
	cp fdelmodc $B

clean:
		rm -f core $(OBJC) $(OBJM) fdelmodc 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/fdelmodc 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



