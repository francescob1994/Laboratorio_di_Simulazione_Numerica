CC = g++
CFLAGS = -Wall -O3 --std=c++11

lsn_01.exe : lsn_01.o random.o
	$(CC) random.o lsn_01.o -o lsn_01.exe
lsn_01.o : lsn_01.cxx 
	$(CC) -c lsn_01.cxx -o lsn_01.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o lsn_01.exe seed.out
