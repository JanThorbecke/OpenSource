# Include general defines and machine dependent stuff

include ../../Make_include

LIBS = -L$L -lextrap3d  -lgenfft -lblas -llapack -lm 

LDFLAGS += $(LIBS)

CPROGS = \
	cfpmod3d \
	extrap3d \
	tablecalc3d \


CTARGET	  = %: %.c

INSTALL: $(CPROGS) 

$(CPROGS):	$(CTARGET) $L/libextrap3d.a
	$(CC) $(CFLAGS) $(OPTC) $(CPPFLAGS) $@.c $(LDFLAGS) -o $@.tmp
	test ! -f $@ || $(RM) $@
	mv $@.tmp $B/$@

clean:
	-rm -f *.o

realclean:
	$(MAKE) clean
	rm -f $(CPROGS) 

