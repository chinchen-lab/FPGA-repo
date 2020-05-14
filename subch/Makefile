#
# Makefile
#
CC=g++
# if you want to use debugger, add -g to CFLAGS and LDFLAGS
CFLAGS=-std=c++17 -O2 -O3 -g
INCLUDES=-I../src/
# INCLUDES=-I../include/
HEADERS=fpga_gr.h node.h
LFLAGS=
# LFLAGS=-L../lib/
#LIBS=-lm -lsystemc
# LIBS=-lm -lgurobi_c++ -lgurobi70 -fopenmp
SOURCES=main.cpp fpga_gr.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../bin/fpga_subch

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJECTS) -o $@ $(LFLAGS) $(LIBS)

%.o:  %.cpp  $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

clean:
	rm -rf *.o $(EXECUTABLE)
