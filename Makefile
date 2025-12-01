CC = gcc
CFLAGS = -g -Wall -Werror -O3
TARGETS = scf_test

.PHONY: default all clean

default: all

all: $(TARGETS)

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ -c $<

scf_test: scf_test.o scf_filter.o scf_tx.o scf_rx.o scf_packet.o
	$(CC) $(CFLAGS) -o $@ $? -lfftw3f -lgsl -lm -lfec

clean:
	-rm -f *.o $(TARGETS)
