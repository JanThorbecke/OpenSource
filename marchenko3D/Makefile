# Makefile

include ../Make_include

LIBS    += -L$L -lgenfft -lzfp -lm $(LIBSM)
CFLAGS  += 
#OPTC += -g -O0 -Wall 

#ALL: fmute marchenko marchenko2

ALL: marchenko3D fmute3D TWtransform

SRCT	= marchenko3D.c \
		getFileInfo3D.c  \
		getFileInfo3DW.c  \
		getFileInfo3Dzfp.c  \
		readData3D.c \
		readShotData3D.c \
		readShotData3Dzfp.c \
		readTinvData3D.c \
		synthesis3D.c \
		applyMute3D.c \
		writeData3D.c \
		makeWindow3D.c \
		ampest3D.c \
		imaging3D.c \
		homogeneousg3D.c \
		readSnapData3D.c \
		writeDataIter3D.c \
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

SRCTW	= TWtransform.c \
		getFileInfo3D.c  \
		writeData3D.c \
		wallclock_time.c \
		name_ext.c  \
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

OBJTW	= $(SRCTW:%.c=%.o)

TWtransform:	$(OBJTW) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o TWtransform $(OBJTW) $(LIBS)

install: marchenko3D fmute3D TWtransform
	cp marchenko3D $B
	cp fmute3D $B
	cp TWtransform $B

clean:
		rm -f core marchenko3D $(OBJT) fmute3D $(OBJJ3) TWtransform $(OBJTW)

realclean: clean
		rm -f $B/marchenko3D $B/fmute3D $B/TWtransform
