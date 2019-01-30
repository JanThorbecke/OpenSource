# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -O0 -Wall 

#ALL: fmute marchenko marchenko2

ALL: fmute marchenko marchenko3D fmute3D

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

SRCT	= marchenko3D.c \
		getFileInfo3D.c  \
		readData3D.c \
		readShotData3D.c \
		readTinvData3D.c \
		synthesis3D.c \
		applyMute.c \
		writeData.c \
		writeDataIter.c \
		wallclock_time.c \
		name_ext.c  \
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

SRCJ3	= fmute3D.c \
		getFileInfo3D.c  \
		readData3D.c \
		applyMute.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

OBJJ	= $(SRCJ:%.c=%.o)

fmute:	$(OBJJ) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute $(OBJJ) $(LIBS)

OBJT	= $(SRCT:%.c=%.o)

marchenko3D:	$(OBJT) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko3D $(OBJT) $(LIBS)

OBJH	= $(SRCH:%.c=%.o)

marchenko:	$(OBJH) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko $(OBJH) $(LIBS)

OBJH2	= $(SRCH2:%.c=%.o)

marchenko2:	$(OBJH2) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko2 $(OBJH2) $(LIBS)

OBJJ3	= $(SRCJ3:%.c=%.o)

fmute3D:	$(OBJJ3) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute3D $(OBJJ3) $(LIBS)

install: fmute marchenko test fmute3D
	cp fmute $B
	cp marchenko $B
	cp marchenko3D $B
	cp fmute3D $B

#	cp marchenko2 $B

clean:
		rm -f core fmute $(OBJJ) marchenko $(OBJH) marchenko2 $(OBJH2) marchenko3D $(OBJT) fmute3D $(OBJJ3)

realclean: clean
		rm -f $B/fmute $B/marchenko $B/marchenko2 $B/marchenko3D $B/fmute3D
