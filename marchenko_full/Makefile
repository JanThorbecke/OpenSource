# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -O0 -Wall 

#ALL: fmute marchenko marchenko2

ALL: fmute marchenko

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
		freqwave.c

OBJJ	= $(SRCJ:%.c=%.o)

fmute:	$(OBJJ) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute $(OBJJ) $(LIBS)

OBJH	= $(SRCH:%.c=%.o)

marchenko:	$(OBJH) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko $(OBJH) $(LIBS)

install: fmute marchenko 
	cp fmute $B
	cp marchenko $B

#	cp marchenko2 $B

clean:
		rm -f core fmute $(OBJJ) marchenko $(OBJH)

realclean: clean
		rm -f $B/fmute $B/marchenko




