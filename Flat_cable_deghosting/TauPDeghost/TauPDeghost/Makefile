CC = g++

BIN = test

SRC= $(wildcard *.C)

OBJS = $(patsubst %.C,%.o,$(SRC))


#FFTWDIR=../../../../lib/fftw

BUILD_OPTS= -g -O3 -fopenmp 

LIB_FLAGS= -lfftw3 -lm 


all: $(BIN) 
	@echo OBJS: $(OBJS)
	
$(BIN): $(OBJS)
	$(CC) $(BUILD_OPTS) $(OBJS) -o $@  -L$(FFTWDIR)/lib  $(LIB_FLAGS)

%.o:%.C
	$(CC) $(BUILD_OPTS) -c $< -o $@    -I$(FFTWDIR)/include

clean:
	rm -rf *.o $(BIN)
