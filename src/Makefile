C++ = g++
CFLAGS = -O2 -Wall -Wextra -Wpedantic -std=gnu++11 -pthread

OBJS = $(BUILD_DIR)/main.o $(BUILD_DIR)/alignment.o $(BUILD_DIR)/ConfigDB.o $(BUILD_DIR)/templateFinder.o $(BUILD_DIR)/birAligner.o $(BUILD_DIR)/birConsolidate.o $(BUILD_DIR)/printOptions.o $(BUILD_DIR)/birFinder.o $(BUILD_DIR)/Fragment.o $(BUILD_DIR)/programExecution.o

BUILD_DIR = Build

m3birfinder: $(OBJS)
	$(C++) -o Bin/m3birfinder $(CFLAGS) $(OBJS) 

$(BUILD_DIR)/main.o: main.cpp defs.h programExecution.h programExecution.cpp
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/main.o -c main.cpp

$(BUILD_DIR)/alignment.o: alignment.cpp defs.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/alignment.o -c alignment.cpp

$(BUILD_DIR)/ConfigDB.o: ConfigDB.cpp ConfigDB.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/ConfigDB.o -c ConfigDB.cpp 

$(BUILD_DIR)/templateFinder.o: templateFinder.cpp defs.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/templateFinder.o -c templateFinder.cpp 

$(BUILD_DIR)/birAligner.o: birAligner.cpp defs.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/birAligner.o -c birAligner.cpp 

$(BUILD_DIR)/birConsolidate.o:
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/birConsolidate.o -c birConsolidate.cpp 

$(BUILD_DIR)/printOptions.o: printOptions.cpp defs.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/printOptions.o -c printOptions.cpp 

$(BUILD_DIR)/birFinder.o: birFinder.cpp
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/birFinder.o -c birFinder.cpp 

$(BUILD_DIR)/Fragment.o: Fragment.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/Fragment.o -c Fragment.cpp 

$(BUILD_DIR)/programExecution.o: defs.h
	$(C++) $(CFLAGS) -o $(BUILD_DIR)/programExecution.o -c programExecution.cpp 

clean:
	rm -rf Bin/*
	rm -rf Build/*

