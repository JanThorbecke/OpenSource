# Makefile for ...FFT

include ../../Make_include

LIB = $L/libextrap2d.a

CFLAGS += -D_GNU_SOURCE 

CSRC = \
	GaussWindow.c \
	KaiserWindow.c \
	getFileInfo.c \
	getModelInfo.c \
	readData.c \
	writeData.c \
	forwExtr.c \
	getrecextr.c \
	getrecvsp.c \
	invExtr.c \
	onewayextr.c \
	kxwfilter.c \
	minmax.c \
	remez.c \
	setzsrc.c \
	shortoper.c \
	toeplitz.c \
	spline3.c \
	srcarray.c \
	tablecalc_opt.c \
	findBestOper.c \
	extrapEdge.c \
	trunc1D.c \
	weightfunct.c \
	xwMigrOpl.c \
	xwZoMigr.c \
	xwBeam.c \
	xwCFP.c \
	xwExtr.c \
	xwSnap.c \
	xwVSP.c \
	xwExtrG.c \
	kwExtr.c \
	kwMigr.c \
	kwZoMigr.c \
	atopkge.c \
	docpkge.c \
	getpars.c \
	verbosepkg.c \
	wallclock_time.c \

#if there is a Fortran compiler you can also use Remez optimised operators
# uncomment this line to do so.
#FSRC = \
#	optRemez.f \

OBJC    = $(CSRC:%.c=%.o)
OBJF    = $(FSRC:%.f=%.o)

CARCH = $(patsubst %.c, $(LIB)(%.o), $(CSRC))
FARCH = $(patsubst %.f, $(LIB)(%.o), $(FSRC))

ARCH := $(CARCH) $(FARCH)

# define all object files
#OBJS = $(patsubst $(LIB)(%), %, $(ARCH))

# give rule to produce library
INSTALL: mkdirs $(LIB)

$(LIB): $(ARCH)
	$(RANLIB) $(LIB)

mkdirs:
	-mkdir $(ROOT)/lib
	-mkdir $(ROOT)/include
	-ln -sf $(ROOT)/extrap/lib/optim.h $I/optim.h 

remake	:
	@-rm -f $(LIB) INSTALL
	@$(MAKE)

list	:
	$(AR) tv $(LIB)

clean:
	-rm -f core a.out *.o $(LIB)

realclean:
	$(MAKE) clean
	-rm -rf $(LIB)

.PRECIOUS:	 $(LIB)

