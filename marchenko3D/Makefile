# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -O0 -Wall 

#ALL: fmute marchenko marchenko2

ALL: marchenko3D fmute3D

SRCT	= marchenko3D.c \
		getFileInfo3D.c  \
		readData3D.c \
		readShotData3D.c \
		readTinvData3D.c \
		synthesis3D.c \
		applyMute3D.c \
		writeData3D.c \
		ampest3D.c \
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
		applyMute3D.c \
		writeData3D.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

OBJT	= $(SRCT:%.c=%.o)

marchenko3D:	$(OBJT) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko3D $(OBJT) $(LIBS)

OBJJ3	= $(SRCJ3:%.c=%.o)

fmute3D:	$(OBJJ3) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fmute3D $(OBJJ3) $(LIBS)

install: marchenko3D fmute3D
	cp marchenko3D $B
	cp fmute3D $B

clean:
		rm -f core marchenko3D $(OBJT) fmute3D $(OBJJ3)

realclean: clean
		rm -f $B/marchenko2 $B/marchenko3D $B/fmute3D
