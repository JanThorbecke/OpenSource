# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC  = -I.
#LIBS    += -L$L -lgenfft -lm $(LIBSM)
#LIBS    += -L$L -lgenfft -lm -lc
#OPTC = -g -Wall -fsignaling-nans -O0
#OPTC += -fopenmp -Waddress
#OPTC += -g -O0
#OPTC := $(subst -O3 -ffast-math, -O1 -g ,$(OPTC))
#PGI options for compiler feedback
#OPTC += -Mprof=lines
#LDFLAGS += -Mprof=lines

#		side.c \
#		corner.c \
#		near_source.c \
#		Grid2Time1.c \

all: raytime3d

PRG = raytime3d

SRCC	= $(PRG).c \
		vidale3d.c \
		src3D.c \
		getParameters3d.c  \
		getWaveletInfo.c  \
		writeSrcRecPos.c  \
		readModel3d.c  \
		getWaveletHeaders.c  \
		verbosepkg.c  \
        getModelInfo3D.c  \
		wallclock_time.c  \
        recvPar3D.c  \
        writesufile.c  \
        name_ext.c  \
		atopkge.c \
		docpkge.c \
		threadAffinity.c \
		getpars.c \
		writeData3D.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) raytime3d.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o raytime3d $(OBJC) $(LIBS)

install: raytime3d
	cp raytime3d $B

clean:
		rm -f core $(OBJC) $(OBJM) raytime3d 

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/raytime3d 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



