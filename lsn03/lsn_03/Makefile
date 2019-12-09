CC = g++
CFLAGS = -Wall -O3 --std=c++11

lsn_03.exe : lsn_03.o random.o
	$(CC) random.o lsn_03.o -o lsn_03.exe
lsn_03.o : lsn_03.cxx
	$(CC) -c lsn_03.cxx -o lsn_03.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o lsn_03.exe seed.out
