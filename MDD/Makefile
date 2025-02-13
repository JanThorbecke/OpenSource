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
LIBS += $(BLAS)

#CFLAGS += -g # -I$(MKLROOT)/include  
#FFLAGS += -g -mcmodel=medium
#LIBS += -lblas -llapack -L$L -lgenfft $(LIBSM) -lc -lm

ifeq ($(FC),ifort)
   LDFLAGS += -nofor-main
else ifeq ($(FC),ifx)
   LDFLAGS += -nofor-main
endif
ifeq ($(FC),gfortran)
   LDFLAGS += -Wno-main
endif

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

SRCF =  lsqrDataModule.F90 \
        zlsqrDataModule.F90 \
        zlsqrblasInterface.F90 \
        zlsqrModule.F90

OBJC	= $(SRCC:%.c=%.o)

OBJF    = $(SRCF:%.F90=%.o)

$(PRG):	$(OBJC) $(OBJF)
	$(FC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o $(PRG) $(OBJC) $(OBJF) $(LIBS) 

install: $(PRG) 
	cp $(PRG) $B

clean:
		rm -f core $(OBJC) $(OBJF) $(OBJM) *.mod $(PRG) 

realclean:
		rm -f core $(OBJC) $(OBJF) $(OBJM) *.mod $(PRG) $B/$(PRG) 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar


