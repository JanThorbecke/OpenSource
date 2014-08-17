# Makefile for ...FFT

include ../../Make_include

LIB = $L/libextrap3d.a

CFLAGS += -D_GNU_SOURCE 

CSRC = \
	xwExtr3d.c \
	xwCFP3d.c \
	getreccfp.c \
	SolveAxb.c \
	W3d.c \
	conv2D.c \
	conv2DMcC.c \
	conv_kxw.c \
	dataFileIO.c \
	forwExtr.c \
	frequency_distribution.c \
	image_condition.c \
	interpolateXY.c \
	shortoper.c \
	tablecalc_1D.c \
	tablecalc_2D.c \
	taperedges.c \
	velocityFileIO.c \
	wlsq2d8c.c \
	wlsq2dc.c \
	write_image.c \
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
INSTALL: $(LIB)

$(LIB): $(ARCH)
	$(RANLIB) $(LIB)

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

