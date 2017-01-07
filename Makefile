#master makefile for OpenSource 

all: mkdirs 
	cd FFTlib		; $(MAKE)
	cd fdelmodc		; $(MAKE) install
	cd utils		; $(MAKE) install
	cd marchenko	; $(MAKE) install

mkdirs:
	-mkdir -p lib
	-mkdir -p include
	-mkdir -p bin

clean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@

realclean:
	cd FFTlib 		; $(MAKE) $@
	cd fdelmodc		; $(MAKE) $@
	cd utils		; $(MAKE) $@
	cd marchenko	; $(MAKE) $@
	rm -f lib/*
	rm -f include/*
	rm -f bin/*
