CC = gcc
CFLAGS = -g -Wall -Werror -O3 -MMD -MP
LDFLAGS =
LDLIBS = -lfftw3f -lgsl -lm

ALSA_CFLAGS := $(shell pkg-config --cflags alsa 2>/dev/null)
ALSA_LIBS := $(shell pkg-config --libs alsa 2>/dev/null)
ifeq ($(ALSA_LIBS),)
ALSA_LIBS = -lasound
endif

COMMON_SRCS = scf_filter.c scf_tx.c scf_rx.c scf_packet.c scf_reed_solomon.c scf_cobs.c
SCF_SRCS = scf_test.c $(COMMON_SRCS)
SCF_ALSA_SRCS = scf_alsa.c $(COMMON_SRCS)
SCF_TEST_ALSA_SRCS = scf_test_alsa.c scf_cobs.c
RS_TEST_SRCS = scf_test_reed_solomon.c

SCF_OBJS = $(SCF_SRCS:.c=.o)
SCF_ALSA_OBJS = $(SCF_ALSA_SRCS:.c=.o)
SCF_TEST_ALSA_OBJS = $(SCF_TEST_ALSA_SRCS:.c=.o)
RS_TEST_OBJS = $(RS_TEST_SRCS:.c=.o)

ALL_OBJS = $(sort $(SCF_OBJS) $(SCF_ALSA_OBJS) $(SCF_TEST_ALSA_OBJS) $(RS_TEST_OBJS))
ALL_DEPS = $(ALL_OBJS:.o=.d)

.PHONY: all clean

all: scf_test scf_test_reed_solomon scf_alsa scf_test_alsa

scf_test: $(SCF_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

scf_test_reed_solomon: scf_test_reed_solomon.o scf_reed_solomon.o
	$(CC) $(LDFLAGS) -o $@ $^

scf_alsa: $(SCF_ALSA_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(ALSA_LIBS) -pthread -lutil

scf_test_alsa: $(SCF_TEST_ALSA_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

scf_alsa.o: CFLAGS += -pthread $(ALSA_CFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

-include $(ALL_DEPS)

clean:
	rm -f scf_test scf_test_reed_solomon scf_alsa scf_test_alsa $(ALL_OBJS) $(ALL_DEPS)
