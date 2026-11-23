CC = gcc
CFLAGS = -g -Wall -Werror -O3 -MMD -MP
LDFLAGS =
LDLIBS = -lfftw3f -lgsl -lm

SCF_SRCS = scf_test.c scf_filter.c scf_tx.c scf_rx.c scf_packet.c scf_reed_solomon.c
RS_TEST_SRCS = scf_test_reed_solomon.c
SCF_OBJS = $(SCF_SRCS:.c=.o)
RS_TEST_OBJS = $(RS_TEST_SRCS:.c=.o)
SCF_DEPS = $(SCF_OBJS:.o=.d)
RS_TEST_DEPS = $(RS_TEST_OBJS:.o=.d)

.PHONY: all clean

all: scf_test scf_test_reed_solomon

scf_test: $(SCF_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

scf_test_reed_solomon: scf_test_reed_solomon.o scf_reed_solomon.o
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

-include $(SCF_DEPS)
-include $(RS_TEST_DEPS)

clean:
	rm -f scf_test scf_test_reed_solomon $(SCF_OBJS) $(RS_TEST_OBJS) $(SCF_DEPS) $(RS_TEST_DEPS)
