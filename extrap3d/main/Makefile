# Include general defines and machine dependent stuff

include ../../Make_include

LIBS += -L$L -lextrap3d -lgenfft -lm 

LDFLAGS += $(LIBS)

CPROGS = \
	cfpmod3d \
	extrap3d \
	tablecalc3d \
	migr3d \
	zomodel3d \

CTARGET	  = %: %.c

INSTALL: $(CPROGS) migr3d_mpi

$(CPROGS):	$(CTARGET) $L/libextrap3d.a
	$(CC) $(CFLAGS) $(OPTC) $(CPPFLAGS) $@.c $(LDFLAGS) -o $@.tmp
	test ! -f $@ || $(RM) $@
	mv $@.tmp $B/$@

migr3d_mpi: migr3d.c  $L/libextrap3d.a
	mpicc -cc=icc $(CFLAGS) -DMPI $(OPTC) $(CPPFLAGS) migr3d.c $(LDFLAGS) -o $@.tmp
	test ! -f $@ || $(RM) $@
	mv $@.tmp $B/$@

clean:
	-rm -f *.o

realclean:
	$(MAKE) clean
	rm -f $(CPROGS) 

