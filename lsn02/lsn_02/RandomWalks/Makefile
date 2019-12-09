CC = g++
CFLAGS = -Wall -O3 --std=c++11

RandomWalks.exe : RandomWalks.o random.o
	$(CC) random.o RandomWalks.o -o RandomWalks.exe
RandomWalks.o : RandomWalks.cxx
	$(CC) -c RandomWalks.cxx -o RandomWalks.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o RandomWalks.exe seed.out
