# Known bugs and issues

Code review findings from 2026-06-10 (Claude Code session). Items 1 and 2 were
confirmed at runtime with ASan/UBSan. Check off items as they get fixed.

## Confirmed crashes / undefined behavior

- [x] **1. Stack buffer overflow in `find_sync` when no signal is present** — `scf_rx.c:129-134`
  *Fixed 2026-06-10: early return from `find_sync` when the first pass finds zero weight; re-verified with ASan.*
  The second refinement pass computes `b0 = max_bin - FFT_RATIO * 2`. If the first
  pass found zero weight everywhere (all-zero demod ring buffer — silence, or the
  first symbols after startup), `max_bin` stays `0`, so `b0 = -8` and the loop
  reads/writes `weight_per_bin[-8]` and `demod_buf[...][negative]`.
  **Verified with ASan**: feeding zero samples crashes with stack-buffer-overflow at
  line 134 on the first `scf_rx()` call.
  Fix: skip the second pass (and threshold check) when the first pass found
  `max_weight == 0`, or clamp `b0` to `cfo_bin_min`.

- [x] **2. NaN converted to `uint8_t` in `demodulate`** — `scf_rx.c:99`
  *Fixed 2026-06-10: `power_sum` is seeded with `FLT_MIN` so the denominator is never
  zero — branchless, keeps the loop vectorizable (an earlier `if (power_sum > 0)`
  guard cost ~35% decode time and was reverted). Re-verified with UBSan
  float-cast-overflow + float-divide-by-zero; decode time back at ~60 ms.*
  `c->demod_buf[demod_buf_idx][i] = WEIGHT_SCALE * power[i] / power_sum;` divides by
  zero when the spectrum is silent, producing NaN; converting NaN to an unsigned
  integer is UB. **Verified with UBSan** (`-fsanitize=float-cast-overflow`):
  `-nan is outside the range of representable values of type 'unsigned char'`.
  Fix: guard against `power_sum == 0` (or add a tiny epsilon).

- [x] **3. `scf.h` does not compile standalone** — `scf.h:29`
  *Fixed 2026-06-10: added `#include <stdbool.h>`; verified a TU including only scf.h compiles with `-Wall -Werror`.*
  The public header uses `bool` in `scf_rx_result` but never includes `<stdbool.h>`.
  It only builds in-repo because `scf_private.h` includes `<stdbool.h>` before
  `scf.h`. Verified: external consumers including `scf.h` directly get
  `unknown type name 'bool'`. Fix: add `#include <stdbool.h>`.

- [x] **4. Always-true assert hides a NULL dereference** — `scf_packet.c:86-97`
  *Fixed 2026-06-10: `packet_type` initialized to `SCF_PACKET_TYPES` and asserted
  `< SCF_PACKET_TYPES` (matching the decode-side assert). Verified the assert fires
  for an unmatched msg_len; scf_test encode path unaffected.*
  `size_t packet_type = -1; ... assert(packet_type >= 0);` — `packet_type` is
  unsigned, the assert can never fire. If `msg_len` matches no packet type,
  `rs_code` stays NULL and `encode_rs_char(NULL, ...)` crashes; also
  `scf_packet_sync_vector[(size_t)-1]` would be indexed.
  Fix: use a signed type or assert `packet_type != (size_t)-1`.

## Logic bugs

- [x] **5. Inverted phase-wrap condition in the TX oscillator** — `scf_tx.c:70-72`
  *Fixed 2026-06-10: condition changed to `< -2.0f * M_PI`. Behavior-neutral as
  expected: full scf_test runs gave 0.49 reception probability before vs 0.54 after
  (n=100, within noise; the suite runs at threshold SNR, ~0.5 is the normal
  operating point).*
  `while (baseband_phase < 2.0f * M_PI) { baseband_phase += 2.0f * M_PI; }`
  should be `< -2.0f * M_PI`. As written, the two wrap loops force the phase into
  [2π, 4π) instead of near zero. Output stays correct only because sinf/cosf are
  periodic, but it defeats the purpose of wrapping and costs float precision.

- [x] **6. `tx_signal` never cleared between test iterations** — `scf_test.c:108`
  *Fixed 2026-06-10: `tx_signal` is now memset alongside `noise` at the top of
  `run_test()`. Full run: 0.59 reception probability (vs 0.49/0.54 before — slight
  improvement consistent with removing ghost-packet interference, though within ~2σ).*
  `run_test()` memsets `noise` but not `tx_signal`. Each test writes its packet at a
  random offset, so the previous iteration's packet remains at its old offset and is
  summed into `rx_signal` as ghost interference, skewing reception probability and
  SNR. Fix: `memset(tx_signal, 0, RX_SIGNAL_LEN * sizeof(float))`.

- [x] **7. `wrong_msg` counter never incremented** — `scf_test.c:176`
  *Fixed 2026-06-10: the "Wrong message" branch now increments `wrong_msg`.
  (Both counters read 0 in a clean run; the counter only registers when RS
  decode succeeds but the payload mismatches.)*
  The "Wrong message" branch increments `wrong_len` (copy-paste error), so the
  `Wrong message:` summary line always prints 0. Fix: `wrong_msg++`.

## Smaller issues

- [x] **8. Re-init doesn't reset state**
  *Fixed 2026-06-10: `scf_rx_init()` now also resets `carrier_phase`,
  `demod_buf_idx`, `input_signal`, and the per-phase demod ring buffers;
  `scf_tx_init()` resets `carrier_phase`, `baseband_phase`, and the modulation
  filter buffer/index. Also added `scf_filter_reset_rx()`/`scf_filter_reset_tx()`
  to clear the overlap-save history buffers in scf_filter.c, called from the
  respective inits. Full suite at 0.58–0.61 reception probability, timing unchanged.*
  `scf_rx_init()` resets `symbol_counter`/`symbol_skip` but not `carrier_phase`,
  `demod_buf_idx`, `input_signal`, or the demod ring buffers; `scf_tx_init()`
  doesn't reset `baseband_phase`, `carrier_phase`, or the modulation filter state.
  Re-initializing to change frequency (explicitly supported per scf.h comment)
  carries stale state; in the test harness RX history leaks across the 100 tests.

- [ ] **9. `scf_filter.c:52`** — the alloc assert checks `fft_out` twice (all five
  buffers happen to be covered anyway). Allocations checked only by `assert`
  vanish under `NDEBUG`.

- [ ] **10. `scf_test.c` printf portability** — `%li`/`%lu` used for `size_t`;
  should be `%zu`. In `dump_signal`, `rv` becomes set-but-unused under `NDEBUG`.

- [ ] **11. `scf.h:1`** — `__SCF_H` is a reserved identifier (double leading
  underscore); prefer `SCF_H`.

## How the runtime bugs were reproduced

Small harness calling `scf_packet_init(42); scf_rx_init(1900.0f);` then `scf_rx()`
with an all-zero `float[SCF_SYM_LEN]` buffer, built with
`-fsanitize=address,undefined` (bug 1) and `-fsanitize=float-cast-overflow` (bug 2).
The shipped `scf_test` never hits these because AWGN keeps every FFT bin nonzero.
