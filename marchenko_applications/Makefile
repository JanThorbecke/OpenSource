# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -Wall

#ALL: fmute marchenko marchenko2

ALL: fmute marchenko_app combine gmshift

SRCJ	= fmute.c \
		getFileInfo.c  \
		readData.c \
		applyMute.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCH	= marchenko.c \
		getFileInfo.c  \
		readData.c \
		readShotData.c \
		readTinvData.c \
		applyMute.c \
		writeData.c \
		writeDataIter.c \
		wallclock_time.c \
		name_ext.c  \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c \
		readSnapData.c \
		Cost.c \
		AmpEst.c \
		freqwave.c \
		getParameters.c \
		getModelInfo.c \
		recvPar.c \
		raytime.c \
		readModel.c \
		JespersRayTracer.c \
		getWaveletHeaders.c \
		getWaveletInfo.c \
		writeSrcRecPos.c \
		writesufile.c \
		gaussGen.c \
		iterations.c \
		imaging.c \
		threadAffinity.c \
		makeWindow.c \
		homogeneousg.c

SRCC	= combine.c \
		getFileInfo.c \
		writeData.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData.c 

SRCG    = gmshift.c \
        getFileInfo.c \
        writeData.c \
        wallclock_time.c \
        getpars.c \
        verbosepkg.c \
        atopkge.c \
        docpkge.c \
        readSnapData.c

OBJJ	= $(SRCJ:%.c=%.o)

fmute:	$(OBJJ) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute $(OBJJ) $(LIBS)

OBJH	= $(SRCH:%.c=%.o)

marchenko_app:	$(OBJH) raytime.h
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko_app $(OBJH) $(LIBS)

OBJC	= $(SRCC:%.c=%.o)

combine:  $(OBJC)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o combine $(OBJC) $(LIBS)

OBJG	= $(SRCG:%.c=%.o)

gmshift:  $(OBJG)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o gmshift $(OBJG) $(LIBS)

install: fmute marchenko_app combine gmshift
	cp fmute $B
	cp marchenko_app $B
	cp combine $B
	cp gmshift $B

#	cp marchenko2 $B

clean:
		rm -f core fmute $(OBJJ) marchenko_app $(OBJH) combine $(OBJC) gmshift $(OBJG)

realclean: clean
		rm -f $B/fmute $B/marchenko_app $B/combine $B/gmshift




