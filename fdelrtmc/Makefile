# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
#LIBS    += -L$L $(LIBSM) -lm -lfftw3 -I/ribarsko/data/holicki/fftw3/fftw-3.3.6-pl1/api -L/ribarsko/data/holicki/fftw3/fftw-3.3.6-pl1/api
LIBS    += -L$L $(LIBSM) -lm 
#LIBS    += -L$L -lgenfft $(LIBSM) -lm 
#LIBS    += -L$L -lgenfft -lm -lc
#OPTC = -g -Wall -fsignaling-nans -O0
#OPTC = -g -O0
#OPTC += -fopenmp -Waddress
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#LDFLAGS += -static

all: fdelrtmc 

PRG = fdelrtmc

SRCC	= $(PRG).c \
		acoustic4.c \
		getParameters.c  \
		fileOpen.c  \
		name_ext.c \
		readModel.c  \
		boundaries.c  \
		verbosepkg.c  \
		wallclock_time.c  \
		atopkge.c \
		docpkge.c \
		getpars.c \
		rtmImagingCondition.c \
		injectSource.c \
		prepareFDOperators.c \
		readSrcWav.c \
		readRcvWav.c \
		readDT.c \
		writeSnapshots.c \
		writesufile.c \
		extractMigrationSnapshots.c \
		createRcvCoordinates.c \
		storeRcvWavefield.c \
		writeMigImage.c \
		writeMigImagePerShot.c \
		writeRec.c \
		decomposition.c \
		fill2DVectorField.c \
		wisdom.c \
		PlaneWaveDecomposition.c \
		k1k2CircFilt.c \
		mvingAvg.c \
		writeModelData.c \
		readModelData.c \
		readRcvCoordinates.c 

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) fdelrtmc.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o fdelrtmc $(OBJC) $(LIBS)

install: fdelrtmc 
	cp fdelrtmc $B

clean:
		rm -f core $(OBJC) $(OBJM) fdelrtmc 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/fdelrtmc 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar

