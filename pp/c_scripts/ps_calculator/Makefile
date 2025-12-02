CC=gcc
CFLAGS=-O2 -Wall -fopenmp
CFLAGS+=-g # Debugging
LIBS=-lm -lfftw3

all: power_spectrum

power_spectrum: power_spectrum.c
	$(CC) $(CFLAGS) -o power_spectrum power_spectrum.c $(LIBS)

clean:
	rm -f power_spectrum

