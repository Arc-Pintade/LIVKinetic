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
	rm -f ./src/*.o bin/* 
scripts/%.o: scripts/%.cpp
	$(CC) -c $(ROOT) -o $@ $< 
%.o: %.cpp
	$(CC) -c $(ROOT) -o $@ $<
ttbar: scripts/ttbar.o $(OBJ)
	$(CC) -o bin/$@ $^ $(ROOT) 
amunu: scripts/amunu.o $(OBJ)
	$(CC) -o bin/$@ $^ $(ROOT)
