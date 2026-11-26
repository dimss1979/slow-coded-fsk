# SCF — Slow Coded FSK modem

An experimental software modem for extremely low-power shortwave communication in
the 13.56 MHz and 27 MHz bands. SCF implements the physical (PHY) layer only: it
turns short byte packets into a very slow, robust FSK signal and back. Its main
counterpart is the **ASFT** tool, which provides the MAC layer on top of this PHY.

The modem works with an audio-rate passband signal (8 kHz sampling, 1900 Hz
carrier by default). An external RF front end translates that signal to the
target band; the receiver side does the inverse conversion.

## Design

The waveform is deliberately slow and redundant to allow simple, low-power
transceivers and very weak-signal reception:

| Parameter | Value |
| --- | --- |
| Modulation | 300-tone MFSK, 5 baud (5 Hz tone spacing) |
| Modem sample rate | 8000 Hz (2000 Hz baseband, 4x decimation) |
| Signal bandwidth | about 1500 Hz around the carrier |
| Inner code | 8 FSK symbols per byte, random codebook, maximum-likelihood decoding |
| Outer code | Reed–Solomon over GF(256), hard decision |
| Framing | 40-symbol pseudorandom sync vector interleaved into every packet |
| Packet types | 8 user bytes in a 16-byte codeword (4 correctable byte errors), or 108 user bytes in a 130-byte codeword (11 correctable byte errors) |
| Useful data rate | about 1.9 bit/s (8-byte packets), 4 bit/s (108-byte packets) |
| TX level | constant envelope, -3 dB RMS full scale |

The receiver downconverts and filters the signal, demodulates it with FFTs at
four symbol timing phases, searches for the sync vector while estimating and
compensating carrier-frequency offset (up to +/-200 Hz), then applies the inner
and outer decoders. The sync vectors, inner codebook and data scrambler are all
derived from a seed, so **both ends must use the same seed**.

## Building

Dependencies: a C compiler, GNU Make, FFTW3 (float), GSL, ALSA development
headers, pthreads and `libutil` (for `openpty`).

```sh
make
```

This produces:

| Program | Purpose |
| --- | --- |
| `scf_alsa` | ALSA sound-card bridge to an asft-compatible pseudo-terminal |
| `scf_test_alsa` | Test client for that pseudo-terminal |
| `scf_test` | Monte-Carlo PHY test over simulated AWGN and a CW interferer |
| `scf_test_reed_solomon` | Reed–Solomon codec self-test |

Use `make clean` to remove build artifacts.

## ALSA bridge

`scf_alsa` connects an ALSA capture/playback device (8 kHz, mono, float samples)
to a pseudo-terminal. Frames arriving on the PTY in the ASFT wire format
(`0x00 | COBS(packet) | 0x00`) are transmitted as SCF packets; decoded packets
are written back to the PTY in the same format. The bridge is packet-length
agnostic — validation belongs to the upper layer. By default it operates in
half duplex and suppresses its own transmission in the captured signal.

```sh
scf_alsa [options]
  -d DEVICE   set both capture and playback device (default: default)
  -C DEVICE   capture device
  -P DEVICE   playback device
  -f HZ       SCF carrier frequency (default: 1900)
  -s SEED     SCF code seed, must match the other end (default: 777)
  -l PATH     PTY symlink path (default: /tmp/scf_alsa)
  -x          full duplex (default: half duplex)
  -h          show help
```

Example:

```sh
scf_alsa -d hw:0 -f 1900 -s 777
```

The companion test client sends one packet or prints received packets:

```sh
scf_test_alsa tx [PTY_PATH]
scf_test_alsa rx [PTY_PATH]
```

`tx` sends a single 8-byte packet (4 random bytes plus CRC32); `rx` receives
frames until interrupted and checks their CRC32. The default PTY path is
`/tmp/scf_alsa`.

## PHY test

`scf_test` runs 100 randomized transmissions with a random start offset,
frequency offset, additive noise and a narrowband interferer, then reports SNR,
Eb/N0, sync events, decoding results and packet reception probability.

```sh
./scf_test
```

## Embedding the modem

The core has no ALSA dependency: `scf_packet_init()`, `scf_tx_init()` and
`scf_rx_init()` set up the codec, then `scf_tx_start()`, `scf_tx()` and
`scf_rx()` operate on one FSK symbol (1600 samples) at a time. See `scf.h` for
the API and `scf_test.c` for a complete send/receive loop.

```c
#include "scf.h"

scf_packet_init(777);          /* seed: must match on both ends */
scf_tx_init(1900.0f);          /* carrier frequency in Hz */
scf_rx_init(1900.0f);

scf_tx_start(msg, 8);          /* only the configured packet lengths work */
while (scf_tx(signal)) {
    /* write signal[0..SCF_SYM_LEN-1] to the audio device */
}

scf_rx(&result, signal);       /* feed one symbol-length block at a time */
if (result.got_msg) {
    /* result.msg / result.msg_len hold the decoded packet */
}
```

Parameters are defined in `scf.h` (tone count, sample rates, symbol length,
decimation, FEC and packet type count). The per-type user and codeword lengths
are in `scf_packet.c`; keep them matched to the packets used by ASFT.

## Status

This is experimental software: the PHY, framing and tests are functional, but
the design is still evolving and is not intended for production use.
