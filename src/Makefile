C++ = g++
CFLAGS = -O2 -Wall -Wextra -Wpedantic

OBJS = main.o alignment.o ConfigDB.o templateFinder.o birAligner.o birConsolidate.o printOptions.o birFinder.o Fragment.o programExecution.o

m3birfinder: $(OBJS)
	$(C++) -o Bin/m3birfinder $(CFLAGS) $(OBJS) 

main.o: main.cpp defs.h
	$(C++) $(CFLAGS) -c main.cpp

alignment.o: alignment.cpp defs.h
	$(C++) $(CFLAGS) -c alignment.cpp

ConfigDB.o: ConfigDB.cpp ConfigDB.h
	$(C++) $(CFLAGS) -c ConfigDB.cpp 

templateFinder.o: templateFinder.cpp defs.h
	$(C++) $(CFLAGS) -c templateFinder.cpp 

birAligner.o: birAligner.cpp defs.h
	$(C++) $(CFLAGS) -c birAligner.cpp 

birConsolidate.o:
	$(C++) $(CFLAGS) -c birConsolidate.cpp 

printOptions.o: printOptions.cpp defs.h
	$(C++) $(CFLAGS) -c printOptions.cpp 

birFinder.o: birFinder.cpp
	$(C++) $(CFLAGS) -c birFinder.cpp 

Fragment.o: Fragment.h
	$(C++) $(CFLAGS) -c Fragment.cpp 

programExecution.o: defs.h
	$(C++) $(CFLAGS) -c programExecution.cpp 

clean:
	rm *.o
	rm -rf Bin/*

