# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
LIBS    += -L$L -lgenfft -lm $(LIBSM)
#LIBS    += -L$L -lgenfft -lm -lc
#OPTC = -g -Wall -fsignaling-nans -O0
#OPTC += -fopenmp -Waddress
OPTC += -g -DDEBUG

PRG = fdemmodc

SRCM	= fdemmodc.c \
		em4.c \
		defineSource.c  \
		getEmParameters.c  \
		getWaveletInfo.c  \
		getModelInfo.c  \
		applySource.c  \
		getRecTimes.c  \
		getBeamTimes.c  \
		writeEmSnapTimes.c  \
		writeEmRec.c  \
		writeSrcRecPos.c  \
		decomposition.c  \
		fileOpen.c  \
		recvPar.c  \
		readEmModel.c  \
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
		getpars.c

OBJM	= $(SRCM:%.c=%.o)

$(PRG): $(OBJM) fdelmodc.h
	$(CC) $(LDFLAGS) $(OPTM) $(CFLAGS) -o fdemmodc $(OBJM) $(LIBS)

install: fdemmodc
	cp fdemmodc $B

clean:
		rm -f core $(OBJC) $(OBJM) fdemmodc

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/fdemmodc


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



