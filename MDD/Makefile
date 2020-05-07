# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
#BLAS libs with Intel compiler
#LIBS += -mkl -L$L -lgenfft $(LIBSM)
#General BLAS library
#LIBS += -L$L -lgenfft $(LIBSM)
#General BLAS library
#LIBS += $(BLAS)

CFLAGS += -I$(MKLROOT)/include
#LIBS += -lblas -llapack -L$L -lgenfft $(LIBSM) -lc -lm

all: mdd 

PRG = mdd

SRCC =  $(PRG).c \
	atopkge.c \
	docpkge.c \
	getpars.c \
	readShotData.c \
	writeEigen.c \
	deconvolve.c \
	computeMatrixInverse.c \
	getFileInfo.c \
	verbosepkg.c \
	name_ext.c \
	wallclock_time.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) 
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o $(PRG) $(OBJC) $(LIBS)

install: $(PRG) 
	cp $(PRG) $B

clean:
		rm -f core $(OBJC) $(OBJM) $(PRG) 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/$(PRG) 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



