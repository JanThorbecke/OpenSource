# Makefile for ...FFT

include ../../Make_include

LIB = libslatec.a
FFLAGS= -O3

FSRC =			\
	d1mach.f \
	dgamln.f \
	i1mach.f \
	zacai.f \
	zacon.f \
	zairy.f \
	zasyi.f \
	zabs.f \
	zbesh.f \
	zbinu.f \
	zbknu.f \
	zbuni.f \
	zbunk.f \
	zdiv.f \
	zexp.f \
	zkscl.f \
	zlog.f \
	zmlri.f \
	zmlt.f \
	zrati.f \
	zs1s2.f \
	zseri.f \
	zshch.f \
	zsqrt.f \
	zuchk.f \
	zunhj.f \
	zuni1.f \
	zuni2.f \
	zunik.f \
	zunk1.f \
	zunk2.f \
	zuoik.f \
	zwrsk.f 
 
OBJF    = $(FSRC:%.f=%.o)

FARCH = $(patsubst %.f, $(LIB)(%.o), $(FSRC))

ARCH := $(FARCH)

install: $(LIB)

$(LIB)  : $(ARCH)
	$(RANLIB) $(LIB)

remake	:
	@-rm -f $(LIB) 
	@$(MAKE)

list	:
	$(AR) tv $(LIB)

clean:
	-rm -f core a.out $(LIB)

realclean:
	$(MAKE) clean
	-rm -f $(LIB) 

.PRECIOUS:	 $(LIB)


