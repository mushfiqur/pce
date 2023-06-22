# this makefile is intended for g++ on Linux

CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -g -Wno-unused-result
LDFLAGS = -lgmpxx -lgmp -lboost_iostreams -lboost_system -lboost_filesystem
# LDFLAGS = -lgmpxx -lgmp

INC = -I/home/mushf/pce/include
SOURCES = $(wildcard ./src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $< -o $@

clean:
	-rm $(OBJECTS)
