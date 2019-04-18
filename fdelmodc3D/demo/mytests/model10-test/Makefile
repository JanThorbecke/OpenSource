#Make cheatsheet
#out.o: src.c src.h
#  $@   # "out.o" (target)
#  $<   # "src.c" (first prerequisite)
#  $^   # "src.c src.h" (all prerequisites)

#gcc flags cheatsheet
#-I($PATH): look for headers here (besides standard dirs). 
#-L($PATH): look for libfiles here (besides standard dirs). Only exe's need libfiles.


MOBJ= checkNaN.o
OBJ= 


all: checkNaN.exe

checkNaN.exe: $(MOBJ) $(OBJ)
	gcc -o $@ $< -L$(CWPROOT)/lib -lsu -lpar -lcwp -lm


%.o: %.c
	gcc -c $< -I$(CWPROOT)/include

clean:
	rm *.o *.exe

