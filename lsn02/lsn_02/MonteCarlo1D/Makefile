CC = g++
CFLAGS = -Wall -O3 --std=c++11

MonteCarlo1D.exe : MonteCarlo1D.o random.o
	$(CC) random.o MonteCarlo1D.o -o MonteCarlo1D.exe
MonteCarlo1D.o : MonteCarlo1D.cxx
	$(CC) -c MonteCarlo1D.cxx -o MonteCarlo1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o MonteCarlo1D.exe seed.out
