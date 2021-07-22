#master makefile for OpenSource 

include Make_include

all: mkdirs 
	cd FFTlib		; $(MAKE)
	cd fdelmodc		; $(MAKE) install
	cd zfp			; $(MAKE) install
	cd utils		; $(MAKE) install
	cd marchenko	; $(MAKE) install
	cd corrvir		; $(MAKE) install
	cd raytime		; $(MAKE) install
ifneq ($(strip $(MKLROOT)),)
	cd fdacrtmc		; $(MAKE) install
else
	@echo "***************************************************************************";
	@echo "**** There is no MKL or other library for the FFTW calls in use by fdacrtmc";
endif
	cd fdelmodc3D	; $(MAKE) install
	cd marchenko3D	; $(MAKE) install
	cd vmar	; $(MAKE) install
ifneq ($(strip $(FC)),)
	cd MDD			; $(MAKE) install
else
	@echo "***************************************************************************";
	@echo "**** There is no Fortran compiler (FC) defined in Make_include to make MDD";
endif

mkdirs:
	-mkdir -p lib
	-mkdir -p include
	-mkdir -p bin

clean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@
	cd corrvir		; $(MAKE) $@
	cd raytime		; $(MAKE) $@
	cd fdacrtmc		; $(MAKE) $@
	cd zfp			; $(MAKE) $@
	cd fdelmodc3D	; $(MAKE) $@
	cd marchenko3D	; $(MAKE) $@
	cd MDD			; $(MAKE) $@
	cd vmar			; $(MAKE) $@

realclean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@
	cd corrvir		; $(MAKE) $@
	cd raytime		; $(MAKE) $@
	cd fdacrtmc		; $(MAKE) $@
	cd zfp			; $(MAKE) $@
	cd fdelmodc3D	; $(MAKE) $@
	cd marchenko3D	; $(MAKE) $@
	cd MDD			; $(MAKE) $@
	cd vmar			; $(MAKE) $@
	rm -f lib/*
	rm -f include/*
	rm -f bin/*
