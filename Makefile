# makefile for poe program
CC=g++
BOOST=/mnt/Data10/rnb203/programs/boost_1_59_0/boost_1_59_0/lib/
INCLUDE=-I/mnt/Data10/rnb203/programs/boost_1_59_0/
FLAGS=-L$(BOOST) -Wl,-rpath,$(BOOST) -lboost_iostreams -std=c++0x -static -O3 -Wall -g
EX_FLAGS=-L/usr/lib/x86_64-redhat-linux5E/lib64/ -lm -L/mnt/Data10/rnb203/programs/anaconda2.2.0/pkgs/zlib-1.2.8-0/lib/ -lz
OB_FLAGS= -c
SOURCES=main.cpp mach.cpp impute2.cpp impute2_haps.cpp general_functions.cpp impute2_trio.cpp impute2_haps_trio.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=poe_generator

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(INCLUDE) $(OBJECTS) $(FLAGS) -o $@ $(EX_FLAGS)

.cpp.o:
	$(CC) $(INCLUDE) $< -o $@ $(FLAGS) $(OB_FLAGS)

clean:
	rm $(OBJECTS) $(EXECUTABLE)
