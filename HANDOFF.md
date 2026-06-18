# Handoff notes — Camera hybrid (TWINS / pump-probe / live-view)

_Last updated 2026-06-15. Folder moved from `C:\Users\mguizzardi\Desktop\Camera hybrid clean` to `D:\Camera hybrid clean` because the C: drive filled to 100%._

## What this app is
PyQt-based control app for a hybrid camera setup doing transmission / differential-transmission
(ΔT) imaging, with several acquisition modes. Entry point is **`main_launcher_lw.py`**, which
opens sub-windows for each mode.

### Key files
- `main_launcher_lw.py` — main launcher window; connects hardware (delay stage, TWINS stage,
  lock-in), opens the sub-mode windows, holds the status LEDs and log.
- `sub_live_lw.py` — live view.
- `sub_deltat_lw.py` — ΔT (differential transmission) static acquisition.
- `sub_kspace_lw.py` — k-space mode.
- `sub_pumpprobe_lw.py` — pump-probe (delay-stage scan).
- `sub_twins_lw.py` — TWINS interferometer static scan.
- `sub_twins_pumpprobe_lw.py` — TWINS + pump-probe combined scan.
- `sub_interval_lw.py`, `sub_lockin_lw.py`, `sub_lockin_pumpprobe_lw.py` — interval / lock-in modes.
- `calibration.py` — shared calibration loaders (motor-jitter correction etc.), used by the
  TWINS static and TWINS pump-probe windows.
- `stage_delay.py`, `stage_driver.py`, `driver_lockin.py`, `labview_manager.py` — hardware drivers.
- `save_config.py` — save-path configuration.

Each frame pair is **odd** (pump off) and **even** (pump on). Derived quantities:
- `Ton`  = `even`              (pump-on transmission)
- `Tavg` = `(odd + even) / 2`
- `DT`   = `even - odd`
- `DT/T` = `(even - odd) / odd * 100`  (zeroed where `|odd| < 1.0` to avoid divide-by-noise)

## Uncommitted changes in this session (NOT yet committed to git)
`git status` will show these 7 files as modified. They have been verified to compile
(`python -m py_compile`). Summary:

### 1. Display/plot modes unified across all acquisition windows
Old combo boxes offered `["DeltaT (dT/T) (%)", "Transmission (T)", "DeltaT (dT)"]`.
Now standardized everywhere to `["DT", "DT/T", "Ton", "Tavg"]`, **default = Tavg (index 3)**.
Index mapping is now: `0=DT`, `1=DT/T`, `2=Ton`, `3=Tavg`.
Affected: `sub_deltat_lw.py`, `sub_kspace_lw.py`, `sub_pumpprobe_lw.py`,
`sub_twins_lw.py`, `sub_twins_pumpprobe_lw.py`.
Background-subtraction logic was updated to match the new indices (DT/DT-T use odd as bg;
Ton/Tavg use the odd/even average as bg).

### 2. Saving now writes ALL four quantities (was: single displayed image)
Saving changed from a single `.npy` of the currently displayed mode to a **`.npz`** bundle
containing `Ton`, `Tavg`, `DT`, `DT_T`, plus `raw_odd` / `raw_even`. The windows now retain the
last odd/even frames (`self._last_odd`, `self._last_even`) so a save can recompute all four.
Filename no longer carries a mode tag (e.g. `<sample>_deltat<pos>_<HHMMSS>.npz`).
There is a fallback that saves only the displayed image if no odd/even pair has been captured yet.
> ⚠️ Consequence: saved files are ~4× larger than before. Make sure the save target has space.

### 3. TWINS connect no longer fakes a connection (`main_launcher_lw.py` `_do_connect_twins`)
Previously, if the real hardware connect failed it set `is_connected = True`, a blue LED, and
logged `[SIM] simulated connection`. This was misleading. Now on failure it:
- sets `is_connected = False`,
- sets the LED red (`#f44336`),
- re-enables the Connect button (retry),
- logs `[FAIL] Twins stage not connected`,
- shows a `QMessageBox.critical` with the underlying exception.

## Known issues / TODO for next agent
- **`_do_connect_delay` in `main_launcher_lw.py` still has the old simulated-connection
  fallback** (blue LED, `[SIM]`). Apply the same real-error treatment as TWINS if desired —
  user was asked, left unchanged for now.
- **Nothing is committed yet.** Once on D: with disk space, run the app to sanity-check, then
  commit (suggested message covering: unified DT/DT-T/Ton/Tavg modes, multi-quantity .npz save,
  TWINS connect error handling).
- Recent committed history (before this session): calibration.py extraction, motor-jitter
  calibration for TWINS static + pump-probe.

## Environment notes
- Windows; primary shell PowerShell, Bash tool also available.
- The C: drive was 100% full (0 bytes) — this corrupted a file mid-write once during the
  session (recovered via `git checkout`). Work from D: going forward; D: had ~650 GB free.
- Verify the `.git` directory came across with the move so history/commit still work.
