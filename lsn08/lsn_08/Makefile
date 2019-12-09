CC = g++
CFLAGS = -Wall -O3 --std=c++11

VMC.exe : VariationalMC.o random.o
	$(CC) random.o VariationalMC.o -o VMC.exe
VariationalMC.o : VariationalMC.cpp
	$(CC) -c VariationalMC.cpp -o VariationalMC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VMC.exe
