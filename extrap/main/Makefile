# Include general defines and machine dependent stuff

include ../../Make_include

LIBS += -L$L -lextrap2d  -lgenfft -lm 

LDFLAGS += $(LIBS)

MPICPROGS = \
	migr_mpi \

CPROGS = \
	migr \
	cfpmod \
	extrap \
	onewvsp \
	opercalc \


CTARGET	  = %: %.c
MPICTARGET= %: %.c

INSTALL: $(CPROGS) $(MPICPROGS)

$(CPROGS):	$(CTARGET) $L/libextrap2d.a
	$(CC) $(CFLAGS) $(OPTC) $(CPPFLAGS) $@.c $(LDFLAGS) -o $@.tmp
	test ! -f $@ || $(RM) $@
	mv $@.tmp $B/$@

$(MPICPROGS): $(MPICTARGET) $L/libextrap2d.a
ifdef MPICC
	@test -d $(@D) || mkdir -p $(patsubst %/,%,$(@D))
	$(MPICC) -DMPI $(CFLAGS) $(OPTC) $(CPPFLAGS) $(@F).c $(LDFLAGS) -o $@.tmp
	test ! -f $@ || $(RM) $@
	mv $@.tmp $B/$@
else
	@echo cannot make $@ because no MPI-C compiler has been defined
endif


clean:
	-rm -f *.o

realclean:
	$(MAKE) clean
	rm -f $(CPROGS) $(MPICPROGS)

