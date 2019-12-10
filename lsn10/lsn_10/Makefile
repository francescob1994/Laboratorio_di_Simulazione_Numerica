CC = g++
CFLAGS = -Wall -O3 --std=c++11

TSP.exe : TSP.o
	$(CC) TSP.o -o TSP.exe
TSP.o : TSP.cpp TSP.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)

clean :
	rm *.o TSP.exe
