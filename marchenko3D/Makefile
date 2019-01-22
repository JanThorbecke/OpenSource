# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -O0 -Wall 

#ALL: fmute marchenko marchenko2

ALL: fmute marchenko test

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

SRCT	= test.c \
		getFileInfo3D.c  \
		readShotData3D.c \
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
		getpars.c

OBJJ	= $(SRCJ:%.c=%.o)

fmute:	$(OBJJ) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute $(OBJJ) $(LIBS)

OBJT	= $(SRCT:%.c=%.o)

test:	$(OBJT) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o test $(OBJT) $(LIBS)

OBJH	= $(SRCH:%.c=%.o)

marchenko:	$(OBJH) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko $(OBJH) $(LIBS)

OBJH2	= $(SRCH2:%.c=%.o)

marchenko2:	$(OBJH2) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko2 $(OBJH2) $(LIBS)

install: fmute marchenko test
	cp fmute $B
	cp marchenko $B
	cp test $B

#	cp marchenko2 $B

clean:
		rm -f core fmute $(OBJJ) marchenko $(OBJH) marchenko2 $(OBJH2) test $(OBJT)

realclean: clean
		rm -f $B/fmute $B/marchenko $B/marchenko2 $B/test
