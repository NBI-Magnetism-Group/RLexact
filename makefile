OBJECTS = RLexact.o RLio.o RLtables.o RLhamil.o RLsymm.o Diagonal.o nr.o RLutil.o  RLcross.o RLlancz.o RLmatrix.o RLsparse.o regc.o
LIBS = -lm 
FLAGS =  -I. -I/usr/include/g++-3 -O3 -gdwarf-2
EXECUTABLE = RLexact
COMPILER = g++
RLexact:	$(OBJECTS)
		$(COMPILER) $(FLAGS) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

RLtables.o:	RLexact.h RLtables.C 
		$(COMPILER) $(FLAGS) -c RLtables.C

RLexact.o:	RLexact.h RLexact.C
		$(COMPILER) $(FLAGS) -c RLexact.C

regc.o:		RLexact.h regc.cpp 
		$(COMPILER) $(FLAGS) -c regc.cpp

RLio.o:		RLexact.h RLio.C
		$(COMPILER) $(FLAGS) -c RLio.C

RLlancz.o:	RLexact.h RLlancz.C
		$(COMPILER) $(FLAGS) -c RLlancz.C

RLmatrix.o:	RLexact.h RLmatrix.C
		$(COMPILER) $(FLAGS) -c RLmatrix.C

RLsymm.o:	RLexact.h RLsymm.C
		$(COMPILER) $(FLAGS) -c RLsymm.C

RLhamil.o:	RLexact.h RLhamil.C
		$(COMPILER) $(FLAGS) -c RLhamil.C

RLutil.o:	RLexact.h RLutil.C
		$(COMPILER) $(FLAGS) -c RLutil.C

Diagonal.o:	RLexact.h Diagonal.C
		$(COMPILER) $(FLAGS) -c Diagonal.C

nr.o:		RLexact.h cnr.h nr.C
		$(COMPILER) $(FLAGS) -c nr.C

RLcross.o:	RLexact.h RLcross.C
		$(COMPILER) $(FLAGS) -c RLcross.C

RLsparse.o:	RLexact.h RLsparse.C
		$(COMPILER) $(FLAGS) -c RLsparse.C

clean:
		rm -f *.o RLexact
