CC = g++

BIN = test

SRC= $(wildcard *.C)

OBJS = $(patsubst %.C,%.o,$(SRC))


FFTWDIR=../../../lib/fftw

BUILD_OPTS= -g -O3 -fopenmp -I$(FFTWDIR)/include -L$(FFTWDIR)/lib  -lfftw3 -lm 


all: $(BIN) 
	@echo OBJS: $(OBJS)
	
$(BIN): $(OBJS)
	g++ $(COMP_FLAG) $(OBJS) -o $@  $(BUILD_OPTS)

%.o:%.C
	g++ -g -O3 -fopenmp -c $< -o $@    -I$(FFTWDIR)/include

clean:
	rm -rf *.o $(BIN)
