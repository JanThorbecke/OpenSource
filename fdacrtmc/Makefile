# Makefile

include ../Make_include

########################################################################
# define general include and system library
CFLAGS  += -I$(ROOT)/zfp/include
LIBS    += -L$L $(LIBSM) -L$(ROOT)/zfp/lib -lzfp

#LIBS    += -L$L -lgenfft $(LIBSM) -lm 
#LIBS    += -L$L -lgenfft -lm -lc
#OPTC = -g -Wall -fsignaling-nans -O0
#OPTC += -g 
#OPTC += -fopenmp -Waddress
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#LDFLAGS += -static

all: fdacrtmc 

PRG = fdacrtmc

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
		fill2DVectorField.c \
		wisdom.c \
		PlaneWaveDecomposition.c \
		k1k2CircFilt.c \
		mvingAvg.c \
		writeModelData.c \
		readModelData.c \
		readRcvCoordinates.c \
		compressMigrationSnapshots.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) fdacrtmc.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o fdacrtmc $(OBJC) $(LIBS)

install: fdacrtmc 
	cp fdacrtmc $B

clean:
		rm -f core $(OBJC) $(OBJM) fdacrtmc 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/fdacrtmc 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar

