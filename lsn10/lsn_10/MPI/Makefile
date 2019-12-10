CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

nTSP.exe : nTSP.o random.o
	$(CC) random.o nTSP.o -o nTSP.exe
nTSP.o : nTSP.cpp nTSP.h
	$(CC) -c nTSP.cpp -o nTSP.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o nTSP.exe seed.out
