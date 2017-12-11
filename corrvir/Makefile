# Makefile for corrvir

include ../Make_include

ALLINC  = -I.
LIBS    += -L$L -lgenfft -lm 
#OPTC += -g

all: corrvir

PRG = corrvir

SRCC	= $(PRG).c \
		getFileInfo.c \
		correlate.c \
		atopkge.c \
		docpkge.c \
		wallclock_time.c  \
		getpars.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) corr.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTC) -o $(PRG)  $(OBJC) $(LIBS)

install: $(PRG) 
	cp $(PRG) $B

clean:
		rm -f core $(OBJC) $(OBJM) $(PRG)

realclean:
		rm -f core $(OBJC) $(OBJM) $(PRG) $B/$(PRG) 


print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



