#master makefile for OpenExtrap

all: mkdirs 
	cd FFTlib		; $(MAKE)
	cd fdelmodc		; $(MAKE) install
	cd fdemmodc		; $(MAKE) install
	cd utils		; $(MAKE) install
	cd marchenko	; $(MAKE) install

mkdirs:
	-mkdir -p lib
	-mkdir -p include
	-mkdir -p bin

clean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd fdemmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@
	find . -name "._*"

realclean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd fdemmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@
	rm -f lib/*
	rm -f include/*
	rm -f bin/*
