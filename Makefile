CFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs --glibs)

BIN=Evt2Root Root2CalMCP Root2Ana Root2SpectraDCs Root2Count

default: $(BIN)

$(BIN): % : %.cpp 
	g++ -o $@ $< -g -Wall -O3 -std=c++11 -lMinuit $(CFLAGS) $(ROOTLIBS)
	
clean:
	rm $(BIN)
