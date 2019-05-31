#master makefile for OpenSource 

all: mkdirs 
	cd FFTlib		; $(MAKE)
	cd fdelmodc		; $(MAKE) install
	cd fdelmodc3D		; $(MAKE) install
	cd utils		; $(MAKE) install
	cd marchenko		; $(MAKE) install
	cd marchenko3D		; $(MAKE) install
	cd corrvir		; $(MAKE) install
	cd raytime		; $(MAKE) install
	cd MDD			; $(MAKE) install

mkdirs:
	-mkdir -p lib
	-mkdir -p include
	-mkdir -p bin

clean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd fdelmodc3D		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko		; $(MAKE) $@
	cd marchenko3D		; $(MAKE) $@
	cd corrvir		; $(MAKE) $@
	cd raytime		; $(MAKE) $@
	cd MDD			; $(MAKE) $@

realclean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd fdelmodc3D		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko		; $(MAKE) $@
	cd marchenko3D		; $(MAKE) $@
	cd corrvir		; $(MAKE) $@
	cd raytime		; $(MAKE) $@
	cd MDD			; $(MAKE) $@
	rm -f lib/*
	rm -f include/*
	rm -f bin/*
