# Makefile

include ../Make_include

#LIBS    += -L$L -lgenfft -lm $(LIBSM)
#OPTC += -g -O0 -Wall 
#OPTC += -g -traceback #-check-pointers=rw -check-pointers-undimensioned


ALL: fmute marchenko marchenko_primaries 

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
		synthesis.c  \
		getFileInfo.c  \
		readData.c \
		readShotData.c \
		readTinvData.c \
		applyMute.c \
		applyMute_tshift.c \
		writeData.c \
		writeDataIter.c \
		wallclock_time.c \
		name_ext.c  \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCP	= marchenko_primaries.c \
		synthesis.c  \
		getFileInfo.c  \
		readData.c \
		readShotData.c \
		readTinvData.c \
		findFirstBreak.c \
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

OBJH	= $(SRCH:%.c=%.o)

marchenko:	$(OBJH) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko $(OBJH) $(LIBS)

OBJP	= $(SRCP:%.c=%.o)

marchenko_primaries:	$(OBJP) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchenko_primaries $(OBJP) $(LIBS)

install: fmute marchenko marchenko_primaries
	cp fmute $B
	cp marchenko $B
	cp marchenko_primaries $B

clean:
		rm -f core fmute $(OBJJ) marchenko $(OBJH) marchenko_primaries $(OBJP)

realclean: clean
		rm -f $B/fmute $B/marchenko $B/marchenko_primaries




