CPLEX_PREFIX=/opt/local

CFLAGS=-D_GNU_SOURCE -std=gnu11 -I$(CPLEX_PREFIX)/include 
CFLAGS_OPT=-O0 -ggdb
#CFLAGS_OPT=-O2

LDFLAGS=-L$(CPLEX_PREFIX)/lib -lm -lcplex


all: lagrange mip

cpphash.o:
	g++ -c cpphash.cpp

lagrange: *.c cpphash.o
	rm -f main.o
	gcc $(CFLAGS) $(CFLAGS_OPT) -I. -DLAGRANGE -c *.c
	g++ $(LDFLAGS) *.o -o lagrange

mip: *.c cpphash.o
	rm -f main.o
	gcc $(CFLAGS) $(CFLAGS_OPT) -I. -DMIP -c *.c
	g++ $(LDFLAGS) *.o -o mip 

clean:
	rm -f mip lagrange
	rm -f *.o
