OBJECTS = RLexact.o RLio.o RLtables.o RLhamil.o RLsymm.o Diagonal.o nr.o RLutil.o  RLcross.o RLlancz.o RLmatrix.o RLsparse.o regc.o
LIBS = -lm 
FLAGS =  -I. -I/usr/include/g++-3 -O3 -gdwarf-2
FLAGS =  -p -I. -I/usr/include/g++-3 -O0 -gdwarf-2 -ansi
FLAGS =  -I. -I/usr/include/g++-3 -O0 -gdwarf-2
OPTFLAGS =  -I. -I/usr/include/g++-3 -O3 -gdwarf-2
#FLAGS =  -p -I. -I/usr/include/g++-3 -O1 -gdwarf-2 -ansi -pedantic
EXECUTABLE = RLexact
COMPILER = g++
RLexact:	$(OBJECTS)
		$(COMPILER) $(FLAGS) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

RLtables.o:	RLexact.h RLtables.C 
		$(COMPILER) $(OPTFLAGS) -c RLtables.C

RLexact.o:	RLexact.h RLexact.C
		$(COMPILER) $(FLAGS) -c RLexact.C

regc.o:		RLexact.h regc.cpp 
		$(COMPILER) $(OPTFLAGS) -c regc.cpp

RLio.o:		RLexact.h RLio.C
		$(COMPILER) $(FLAGS) -c RLio.C

RLlancz.o:	RLexact.h RLlancz.C
		$(COMPILER) $(OPTFLAGS) -c RLlancz.C

RLmatrix.o:	RLexact.h RLmatrix.C
		$(COMPILER) $(OPTFLAGS) -c RLmatrix.C

RLsymm.o:	RLexact.h RLsymm.C
		$(COMPILER) $(OPTFLAGS) -c RLsymm.C

RLhamil.o:	RLexact.h RLhamil.C
		$(COMPILER) $(OPTFLAGS) -c RLhamil.C

RLutil.o:	RLexact.h RLutil.C
		$(COMPILER) $(OPTFLAGS) -c RLutil.C

Diagonal.o:	RLexact.h Diagonal.C
		$(COMPILER) $(OPTFLAGS) -c Diagonal.C

nr.o:		RLexact.h cnr.h nr.C
		$(COMPILER) $(OPTFLAGS) -c nr.C

RLcross.o:	RLexact.h RLcross.C
		$(COMPILER) $(OPTFLAGS) -c RLcross.C

RLsparse.o:	RLexact.h RLsparse.C
		$(COMPILER) $(OPTFLAGS) -c RLsparse.C

clean:
		rm -f *.o RLexact
