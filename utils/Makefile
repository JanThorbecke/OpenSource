# Makefile

include ../Make_include

#LIBS    += -L$L -lgenfft -lm $(LIBSM)
LIBS    += -L$L -lgenfft -lzfp
#OPTC += -openmp 
#OPTC += -g -O0

ALL: makemod makewave extendModel fconv correigen green green3D basop syn2d mat2su ftr1d mutesnap padmodel truncate combine combine_induced reshape_su HomG snap2shot makeR1D pwshift convhomg

SRCM	= \
		makemod.c  \
		wallclock_time.c  \
		diffraction.c \
		randdf.c \
		elipse.c \
		fractint.c \
		grid.c \
		gridabove.c \
		interpolation.c \
		linearint.c \
		plotexample.c \
		polint.c \
		roughint.c \
		sinusint.c \
		spline.c \
		name_ext.c  \
		allocs.c  \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCW	= \
		makewave.c \
		wallclock_time.c \
		freqwave.c \
		allocs.c  \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCE	= \
		extendModel.c \
		wallclock_time.c \
		allocs.c  \
		getModelInfo.c  \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCF	= fconv.c \
		getFileInfo.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCC	= correigen.c \
		getFileInfo.c  \
		readData.c \
		writeData.c \
        name_ext.c  \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c


SRCG	= green.c \
		getFileInfo.c  \
		getrecpos.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCG3	= green3D.c \
		getFileInfo.c  \
		getrecpos3D.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCB	= basop.c \
		getFileInfo.c  \
		kxwfilter.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCS	= syn2d.c \
		getFileInfo.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCA	= mat2su.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCT	= ftr1d.c \
		getFileInfo.c  \
		readData.c \
		writeData.c \
		wallclock_time.c \
		verbosepkg.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

SRCMS   = mutesnap.c \
        getFileInfo3D.c \
        writeData3D.c \
        verbosepkg.c  \
        getpars.c \
        wallclock_time.c \
        atopkge.c \
        docpkge.c \
        readSnapData3D.c

SRCPM	= padmodel.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCTR	= truncate.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCCO	= combine.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCMR	= makeR1D.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCCI	= combine_induced.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCRS   = reshape_su.c \
        getFileInfo3D.c \
        writeData3D.c \
        getpars.c \
        verbosepkg.c \
        atopkge.c \
        docpkge.c \
        readSnapData3D.c

SRCHG   = HomG.c \
        getFileInfo3D.c \
        readData.c \
        writeData3D.c \
        verbosepkg.c  \
        getpars.c \
        wallclock_time.c \
        atopkge.c \
        docpkge.c \
        readSnapData3D.c

SRCSS   = snap2shot.c \
        getFileInfo3D.c \
        writeData3D.c \
        verbosepkg.c  \
        getpars.c \
        wallclock_time.c \
        atopkge.c \
        docpkge.c \
		name_ext.c \
        readSnapData3D.c

SRCPW	= pwshift.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

SRCCH	= convhomg.c \
		getFileInfo3D.c \
		writeData3D.c \
		wallclock_time.c \
		getpars.c \
		verbosepkg.c \
		atopkge.c \
        docpkge.c \
		readSnapData3D.c 

OBJM	= $(SRCM:%.c=%.o)

makemod:	$(OBJM) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o makemod $(OBJM) $(LIBS)

OBJW	= $(SRCW:%.c=%.o)

makewave:	$(OBJW) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o makewave $(OBJW) $(LIBS)

OBJE	= $(SRCE:%.c=%.o)

extendModel:	$(OBJE) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o extendModel $(OBJE) $(LIBS)

OBJF	= $(SRCF:%.c=%.o)

fconv:	$(OBJF) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o fconv $(OBJF) $(LIBS)

OBJC	= $(SRCC:%.c=%.o)

correigen:	$(OBJC) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o correigen $(OBJC) $(LIBS)

OBJG	= $(SRCG:%.c=%.o)

green:	$(OBJG) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o green $(OBJG) $(LIBS)

OBJG3	= $(SRCG3:%.c=%.o)

green3D:	$(OBJG3) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o green3D $(OBJG3) $(LIBS)

OBJB	= $(SRCB:%.c=%.o)

basop:	$(OBJB) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o basop $(OBJB) $(LIBS)

OBJS	= $(SRCS:%.c=%.o)

syn2d:	$(OBJS) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o syn2d $(OBJS) $(LIBS)

OBJA	= $(SRCA:%.c=%.o)

mat2su:	$(OBJA) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o mat2su $(OBJA) $(LIBS)

OBJT	= $(SRCT:%.c=%.o)

ftr1d:	$(OBJT) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o ftr1d $(OBJT) $(LIBS)

OBJMS	= $(SRCMS:%.c=%.o)

mutesnap:    $(OBJMS)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o mutesnap $(OBJMS) $(LIBS)

OBJPM	= $(SRCPM:%.c=%.o)

padmodel:  $(OBJPM)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o padmodel $(OBJPM) $(LIBS)

OBJTR	= $(SRCTR:%.c=%.o)

truncate:  $(OBJTR)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o truncate $(OBJTR) $(LIBS)

OBJCO	= $(SRCCO:%.c=%.o)

combine:  $(OBJCO)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o combine $(OBJCO) $(LIBS)

OBJMR	= $(SRCMR:%.c=%.o)

makeR1D:  $(OBJMR)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o makeR1D $(OBJMR) $(LIBS)

OBJCI	= $(SRCCI:%.c=%.o)

combine_induced:  $(OBJCI)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o combine_induced $(OBJCI) $(LIBS)

OBJRS    = $(SRCRS:%.c=%.o)

reshape_su:  $(OBJRS)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o reshape_su $(OBJRS) $(LIBS)

OBJHG   = $(SRCHG:%.c=%.o)

HomG:   $(OBJHG)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o HomG $(OBJHG) $(LIBS)

OBJSS   = $(SRCSS:%.c=%.o)

snap2shot:   $(OBJSS)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o snap2shot $(OBJSS) $(LIBS)

OBJPW	= $(SRCPW:%.c=%.o)

pwshift:  $(OBJPW)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o pwshift $(OBJPW) $(LIBS)

OBJCH	= $(SRCCH:%.c=%.o)

convhomg:  $(OBJCH)
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o convhomg $(OBJCH) $(LIBS)

install: makemod makewave extendModel fconv correigen green green3D basop syn2d mat2su ftr1d mutesnap padmodel truncate combine combine_induced reshape_su HomG snap2shot makeR1D pwshift convhomg
	cp makemod $B
	cp makewave $B
	cp extendModel $B
	cp fconv $B
	cp correigen $B
	cp green $B
	cp green3D $B
	cp basop $B
	cp syn2d $B
	cp mat2su $B
	cp ftr1d $B
	cp mutesnap $B
	cp padmodel $B
	cp truncate $B
	cp combine $B
	cp combine_induced $B
	cp reshape_su $B
	cp HomG $B
	cp snap2shot $B
	cp makeR1D $B
	cp pwshift $B
	cp convhomg $B

clean:
		rm -f core $(OBJM) makemod $(OBJW) makewave $(OBJE) extendModel $(OBJF) fconv $(OBJG) $(OBJC) correigen green $(OBJG3) green3D $(OBJB) basop $(OBJJ) syn2d $(OBJS) mat2su $(OBJA) ftr1d $(OBJT) mutesnap $(OBJMS) truncate $(OBJTR) padmodel $(OBJPM) combine $(OBJCO) makeR1D $(OBJMR) reshape_su $(OBJRS) combine_induced $(OBJCI) HomG $(OBJHG) snap2shot $(OBJSS) pwshift $(OBJPW) convhomg $(OBJCH)

realclean: clean
		rm -f $B/makemod $B/makewave $B/extendModel $B/fconv $B/correigen $B/green $B/green3D $B/basop $B/syn2d $B/mat2su $B/ftr1d $B/mutesnap $B/padmodel $B/truncate $B/combine $B/combine_induced $B/reshape_su $B/HomG $B/snap2shot $B/makeR1D $B/pwshift $B/convhomg
