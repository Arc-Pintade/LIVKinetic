CC = g++
ROOT = `root-config --cflags --glibs --ldflags`
SRC=$(shell ls ./src/*.cpp)
OBJ=$(SRC:.cpp=.o)
.PHONY: clean
.PHONY: clean-all
all: ttbar
clean:
	rm -f ./src/*.o
clean-all:
	rm -f ./src/*.o ttbar
%.o: %.cpp
	$(CC) -c $(ROOT) -o $@ $<
ttbar: $(OBJ) 
	$(CC) $^ $(ROOT) -o $@ 
