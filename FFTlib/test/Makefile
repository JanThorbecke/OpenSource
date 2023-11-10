# Makefile

include ../../Make_include

LIBS += -L$L -lgenfft $(FFT) -lm

ALL: cc1test cc2dtest ccmtest rc1test rc2dtest rcmtest rc1loop cc1loop mkl3d

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

rc1loop: rc1loop.o
	$(CC) $(OPTC) $(CFLAGS) -o rc1loop rc1loop.o $(LIBS)

cc1loop: cc1loop.o
	$(CC) $(OPTC) $(CFLAGS) -o cc1loop cc1loop.o $(LIBS)

mkl3d: mkl3d.o
	$(CC) $(OPTC) $(CFLAGS) -o mkl3d mkl3d.o $(LIBS)

clean:
		rm -f core *.o cc1test cc2dtest ccmtest rc1test rc2dtest rcmtest rc1loop cc1loop mkl3d

realclean: clean




