# Makefile for ...FFT

FC = gfortran-mp-13
LIB = libslatec.a
FFLAGS= -ffree-form -O2
FFLAGS= -O2

RANLIB = ar -s
#LIB = $L/libdgenfft.a
#OPTC += -DDOUBLE

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
 
#	xgetua.f \
#	fdump.f \
#	xerhlt.f \
#	xercnt.f \
#	xersve.f \
#	xermsg.f \
#	xerprn.f \
#	j4save.f \

OBJF    = $(FSRC:%.f=%.o)

FARCH = $(patsubst %.f, $(LIB)(%.o), $(FSRC))

ARCH := $(FARCH)

INSTALL: $(LIB)

$(LIB)  : $(ARCH)
	$(RANLIB) $(LIB)

remake	:
	@-rm -f $(LIB) INSTALL
	@$(MAKE)

list	:
	$(AR) tv $(LIB)

clean:
	-rm -f core a.out $(LIB)

#	@cd test ; $(MAKE) clean

realclean:
	$(MAKE) clean
	-rm -f $(LIB) 
#	@cd test ; $(MAKE) realclean 

.o.a    :
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o
.f.o    :
	$(FC) -c $(FFLAGS) $(OPTF) $<

.PRECIOUS:	 $(LIB)


