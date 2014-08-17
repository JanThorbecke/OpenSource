# Makefile

include ../../Make_include

LIBS    += -L$L -lgenfft -lm

ALL: cc1test cc2dtest ccmtest rc1test rc2dtest rcmtest 

cc1test: cc1test.o
	$(CC) $(OPTC) $(CFLAGS) -o cc1test cc1test.o $(LIBS)

cc2dtest: cc2dtest.o
	$(CC) $(OPTC) $(CFLAGS) -o cc2dtest cc2dtest.o $(LIBS)

ccmtest: ccmtest.o
	$(CC) $(OPTC) $(CFLAGS) -o ccmtest ccmtest.o $(LIBS)

rc1test: rc1test.o
	$(CC) $(OPTC) $(CFLAGS) -o rc1test rc1test.o $(LIBS)

rc2dtest: rc2dtest.o
	$(CC) $(OPTC) $(CFLAGS) -o rc2dtest rc2dtest.o $(LIBS)

rcmtest: rcmtest.o
	$(CC) $(OPTC) $(CFLAGS) -o rcmtest rcmtest.o $(LIBS)

clean:
		rm -f core *.o cc1test cc2dtest ccmtest rc1test rc2dtest rcmtest 

realclean: clean




