CPLEX_INCL=/usr/local/include/ilcplex
CPLEX_LIB=/usr/local/lib/libcplex.a 
#/usr/local/lib/libcplexdistmip.a -ldl
#CPLEX_INCL=/opt/local/include
#CPLEX_LIB=/opt/local/lib

#CFLAGS=-fopenmp=libiomp5 -D_GNU_SOURCE -std=gnu11 -I$(CPLEX_INC)
CFLAGS=-D_GNU_SOURCE -std=gnu11 -I$(CPLEX_INC)
#CFLAGS_OPT=-O0 -ggdb
CFLAGS_OPT= -ffast-math -march=native -Ofast

#LDFLAGS=$(CPLEX_LIB) -lm -fopenmp=libiomp5
#LDFLAGS=$(CPLEX_LIB) -lm
LDFLAGS=$(CPLEX_LIB) -ffast-math -Ofast -march=native -lm -lglpk


all: lagrange mip findf gr lp mipfm

cpphash.o:
	clang++ $(CFLAGS_OPT) -c cpphash.cpp

findf: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DLAGRANGE_FINDF -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o findf

lagrange: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DLAGRANGE -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o lagrange

gr: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DGR -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o gr

lp: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DLP -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o lp

lpi: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DLPI -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o lpi

mip: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DMIP -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o mip

mipfm: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) $(CFLAGS_OPT) -I. -DMIPFM -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o mipfm

mipfm2: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) -O0 -ggdb -I. -DMIPFM -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o mipfm2

mipfm3: *.c cpphash.o
	rm -f main.o
	clang $(CFLAGS) -O0 -ggdb -I. -DMIPFM -c *.c
	clang++ *.o $(LDFLAGS) -lpthread -o mipfm3

clean:
	rm -f mip lagrange
	rm -f *.o
