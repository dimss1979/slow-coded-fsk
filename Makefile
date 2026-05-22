CC = gcc
CFLAGS = -g -Wall -Werror -O3 -MMD -MP
LDFLAGS =
LDLIBS = -lfftw3f -lgsl -lm -lfec

SCF_SRCS = scf_test.c scf_filter.c scf_tx.c scf_rx.c scf_packet.c
SCF_OBJS = $(SCF_SRCS:.c=.o)
SCF_DEPS = $(SCF_OBJS:.o=.d)

.PHONY: all clean

all: scf_test

scf_test: $(SCF_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

-include $(SCF_DEPS)

clean:
	rm -f scf_test $(SCF_OBJS) $(SCF_DEPS)
