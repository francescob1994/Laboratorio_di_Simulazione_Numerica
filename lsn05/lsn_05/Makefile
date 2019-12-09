CC = g++
CFLAGS = -Wall -O3 --std=c++11

lsn_05.exe : random.o funzioni.o lsn_05.o
	$(CC) random.o funzioni.o lsn_05.o -o lsn_05.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

funzioni.o : funzioni.cpp funzioni.h random.h 
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

lsn_05.o : lsn_05.cpp funzioni.h random.h 
	$(CC) -c lsn_05.cpp -o lsn_05.o $(CFLAGS)

clean :
	rm *.o lsn_05.exe *.dat
