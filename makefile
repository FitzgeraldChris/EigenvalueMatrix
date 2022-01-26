OPTS = -Wall -std=c++14

all: MatrixMaker

MatrixMaker: Matrix.o main.o
	g++ Matrix.o main.o -o MatrixMaker

main.o: main.cpp Matrix.hpp
	g++ $(OPTS) -c main.cpp

Matrix.o: Matrix.cpp Matrix.hpp
	g++ $(OPTS) -c Matrix.cpp
	
clean:
	rm *.o MatrixMaker