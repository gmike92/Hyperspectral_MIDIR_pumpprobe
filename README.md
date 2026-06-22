# Hyperspectral MID-IR Pump–Probe — Camera Hybrid Control Suite

A Python + LabVIEW control suite for time-resolved mid-infrared imaging and
spectroscopy. A 2D-MCT focal-plane camera is read out by LabVIEW; Python (PyQt6 /
pyqtgraph) drives the experiment: it moves the delay and interferometer stages,
streams and processes frames, and saves the data.

It supports differential-transmission (ΔT / ΔT·T⁻¹) imaging, delay-stage
pump–probe scans, and **TWINS** Fourier-transform hyperspectral acquisition
(static and pump–probe), all sharing one camera and one region of interest.

---

## Table of contents
1. [Concept and architecture](#1-concept-and-architecture)
2. [Hardware](#2-hardware)
3. [Installation](#3-installation)
4. [Repository layout](#4-repository-layout)
5. [Running the app](#5-running-the-app)
6. [The acquisition model: odd/even frames and derived quantities](#6-the-acquisition-model)
7. [The LabVIEW bridge](#7-the-labview-bridge)
8. [The shared ROI](#8-the-shared-roi)
9. [Measurement windows](#9-measurement-windows)
10. [TWINS spectral processing (the pipeline in detail)](#10-twins-spectral-processing)
11. [Data and file formats](#11-data-and-file-formats)
12. [Configuration and persisted settings](#12-configuration-and-persisted-settings)
13. [Troubleshooting](#13-troubleshooting)
14. [Developer notes](#14-developer-notes)

---

## 1. Concept and architecture

The suite uses a **"puppeteer"** split:

- **Python** owns the GUI, scan logic, plotting, signal processing, and file I/O.
- **LabVIEW** owns the low-level camera: a single VI, `Experiment_manager.vi`,
  runs continuously and acquires frames on command.

Python drives the VI over **ActiveX/COM** (`win32com`, from `pywin32`). It writes
a command enum, waits for the VI to return to *Idle*, and reads back the acquired
frames. See [§7](#7-the-labview-bridge).

```
                         ┌─────────────────────────────────────────────┐
                         │  main_launcher_lw.py  (control panel)        │
                         │  • Start Manager / Init Camera               │
                         │  • Connect Delay / TWINS stages              │
                         │  • Open measurement windows                  │
                         └───────────────┬─────────────────────────────┘
                                         │ opens
          ┌──────────────┬───────────────┼───────────────┬──────────────┐
          ▼              ▼               ▼               ▼              ▼
      Live View      ΔT + Stage     Pump–Probe       TWINS FTIR    TWINS Pump–Probe
    sub_live_lw     sub_deltat_lw  sub_pumpprobe_lw  sub_twins_lw  sub_twins_pumpprobe_lw
          │              │               │               │              │
          └──────────────┴───────┬───────┴───────────────┴──────────────┘
                                 │ all share
                  ┌──────────────┴───────────────┐
                  ▼                               ▼
          ROIState (roi_state.py)        LabVIEWManager (labview_manager.py)
          one ROI for every window       ActiveX → Experiment_manager.vi → camera
```

The camera frame stream and the derived quantities are described in
[§6](#6-the-acquisition-model). The optical-delay and interferometer stages are
controlled by dedicated drivers ([§2](#2-hardware)).

---

## 2. Hardware

| Subsystem | Device | Driver / interface | Module |
|---|---|---|---|
| Camera | 2D-MCT focal-plane array | LabVIEW VIs + vendor SDK (`Camera/pt_2dmct1_*.dll`) | `labview_manager.py` |
| Optical delay | Thorlabs Kinesis translation stage | `pylablib` (Kinesis) | `stage_delay.py` |
| Interferometer | **NIREOS GEMINI** (SmarAct) birefringent delay | `SCU3DControl.dll` via `ctypes` | `stage_driver.py` |
| Lock-in | Stanford Research **SR865A** | `qcodes` (VISA/USB) | `driver_lockin.py` |

**Delay convention** (double-pass): `Δt [fs] = 2·Δx [mm] / c`, with
`c = 0.000299792458 mm/fs`. Used everywhere fs↔mm appears.

Any stage can be absent — the launcher disables the corresponding controls and
the camera-only windows still work. The delay stage has a *simulated* fallback
(blue LED) for UI testing without hardware; the TWINS stage does **not** fake a
connection (failed connect ⇒ red LED + error dialog).

---

## 3. Installation

### 3.1 Python
64-bit **Python 3.8+** (developed on 3.12). Core packages:

```bash
pip install -r requirements.txt          # PyQt6 pyqtgraph numpy pandas scipy
```

Plus the hardware/bridge packages, by what you actually use:

```bash
pip install pywin32      # REQUIRED for the LabVIEW ActiveX bridge
pip install pylablib     # delay stage (Thorlabs Kinesis)
pip install qcodes       # lock-in modes (SR865A)
```

> `requirements.txt` lists only the pure-Python core. `pywin32` is required for
> any camera acquisition; `pylablib`/`qcodes` are optional per mode. If a driver
> import fails the launcher prints a `[WARN]` and disables that feature rather
> than crashing.

### 3.2 LabVIEW + camera SDK
- A LabVIEW runtime/IDE able to open `Experiment_manager.vi`.
- The camera vendor SDK in `Camera/` (`pt_2dmct1_32.dll` / `_64.dll` + `.lib`/`.h`)
  must be reachable by the VIs.
- The NIREOS `SCU3DControl.dll` must be on the system path for the TWINS stage.

> **After moving the project folder**, LabVIEW may not find subVIs (it stores
> subVI links by absolute path). See [§13](#13-troubleshooting).

### 3.3 Calibration files (required for TWINS spectra)
Two tab-separated files under `Twins/ASRC calibration/` (tracked in this repo):

- `parameters_cal.txt` — spectral calibration: wavelength (µm) ↔ stage reciprocal
  (1/mm). Maps the stage-Fourier frequency axis onto real wavelengths.
- `parameters_int.txt` — a reference interferogram of a known source used to
  remove the stage motor's reproducible nonlinearity (see [§10.2](#102-motor-jitter-position-calibration)).

If missing, the app warns and falls back (simple `1/freq` mapping, no motor
correction). Loaders live in `calibration.py` and cache once per process.

---

## 4. Repository layout

| Path | Purpose |
|---|---|
| **`main_launcher_lw.py`** | Entry point. Control panel: start LabVIEW, init camera, connect stages, open windows, set save directory. |
| `labview_manager.py` | `LabVIEWManager` singleton — ActiveX bridge to `Experiment_manager.vi`. |
| `stage_delay.py` | `DelayStageDriver` — Thorlabs delay stage (move/home, mm↔fs). |
| `stage_driver.py` | `StageDriver` — NIREOS GEMINI interferometer stage. |
| `driver_lockin.py` | `LockInDriver` — SR865A lock-in (qcodes). |
| `calibration.py` | Shared spectral + motor-jitter calibration loaders. |
| `roi_state.py` | `ROIState` singleton — the one shared ROI, persisted to disk. |
| `roi_readout.py` | Small auto-refreshing ROI-bounds label for every window. |
| `save_config.py` | `SaveConfig` singleton — base save directory + `YYYY/MM/DD` folders. |
| **Measurement windows** | |
| `sub_live_lw.py` | Live View (continuous stream, ROI/pixel analysis). |
| `sub_deltat_lw.py` | ΔT live view + manual delay-stage control. |
| `sub_pumpprobe_lw.py` | Delay-stage pump–probe scan. |
| `sub_twins_lw.py` | TWINS FTIR static hyperspectral scan (+ `SpectrumProcessor`). |
| `sub_twins_pumpprobe_lw.py` | TWINS **pump–probe** hyperspectral scan (own `SpectrumProcessor` with phase correction). |
| `sub_kspace_lw.py` | K-space: keep **every** ROI pixel → per-pixel spectra. |
| `sub_lockin_pumpprobe_lw.py` | Pump–probe via lock-in (no camera). |
| `sub_lockin_lw.py`, `sub_interval_lw.py` | Helper widgets (lock-in control; 3-zone delay intervals). |
| **`Camera/`** | LabVIEW acquisition/ΔT VIs + 2D-MCT camera SDK (`pt_2dmct1_*`). |
| `Experiment_manager.vi` | The persistent manager VI Python drives. |
| `Frame_acquisition_server.vi`, `cameraclose.vi`, `Control 2.ctl` | Supporting VIs / typedef. |
| `Twins/ASRC calibration/` | Calibration data (`parameters_cal.txt`, `parameters_int.txt`, Tamosauskas n(λ)). |

Not tracked (gitignored): scan data (`*.npz/.npy/.csv`), `live_roi_state.json`
(runtime ROI), Qt `*.ini` settings, `__pycache__/`.

---

## 5. Running the app

```bash
python main_launcher_lw.py
```

In the launcher:

1. **🚀 Start Manager** — connect to (or launch) LabVIEW and run
   `Experiment_manager.vi`. LED → green on success.
2. **📷 Init Camera** — open the sensor inside the VI (enum `Init`). LED → green.
3. **Connect stages** — *Connect Delay Stage* (then *Home Delay*), and/or
   *Connect Twins Stage*, as needed by your experiment.
4. **Set the save directory** (defaults to `D:\pumpprobedata`).
5. **Open a measurement window** from *Sub-Windows* ([§9](#9-measurement-windows)).

Closing the launcher closes all sub-windows, disconnects stages, and shuts the
manager down cleanly (camera *Close* → VI stop).

---

## 6. The acquisition model

The camera returns frames in pairs synchronised to the pump modulation:

- **odd** = pump **off**
- **even** = pump **on**

Every window derives the same four quantities from an `(odd, even)` pair:

| Quantity | Definition | Meaning |
|---|---|---|
| `Ton` | `even` | pump-on transmission |
| `Tavg` | `(odd + even) / 2` | average transmission |
| `DT` | `even − odd` | differential transmission |
| `DT/T` | `(even − odd) / odd × 100` | relative ΔT/T (%), zeroed where `|odd| < 1` to avoid divide-by-noise |

The display selector in every window is `["DT", "DT/T", "Ton", "Tavg"]`
(index `0..3`), **default `Tavg`**. An optional global background (a blocked-beam
"dark" pair, stored as separate odd/even frames) is subtracted before the
quantities are formed.

---

## 7. The LabVIEW bridge

`LabVIEWManager` (singleton) talks to `Experiment_manager.vi`, which runs a
`While`-loop + `Case` structure. Python sets the **`Enum`** control to a command
and polls until the VI sets it back to `Idle`.

| Enum | Command | Effect |
|---|---|---|
| 0 | `Idle` | wait |
| 1 | `Init` | open camera, return to Idle |
| 2 | `Getframe` | **continuous** acquire loop (Live View); exit via `Stoplive` |
| 3 | `Close` | close camera |
| 4 | `Measure` | acquire **once** (N frames) → Idle — used by all scans |

**Controls** written by Python: `N` (frames to average), `Acq Trigger`,
`Stoplive`, `end` (stop the manager VI).
**Indicators** read by Python: `Odd`, `Even` (frame pair), `T`.

Scans use `Measure` (single-shot) so Python can read the result deterministically;
Live View uses `Getframe` and just peeks at `Odd`/`Even` each tick.

---

## 8. The shared ROI

The region used to reduce images to a signal is defined **once in Live View** and
reused by every scan window, so the measured intensity is consistent across modes.

- `roi_state.py` — `ROIState` singleton holds the ROI rectangle, the ROI/Pixel
  toggle, and the selected pixel. It is **persisted to `live_roi_state.json`**, so
  the selection survives closing/reopening Live View *and* restarting the app.
- Live View **writes** the store whenever you move/resize the ROI, toggle mode, or
  click a pixel. Every scan window **reads** `get_roi_bounds()` / `sel_row` /
  `sel_col` from it — independent of window open order.
- `roi_readout.py` adds a small auto-refreshing label (`Shared ROI rows r0:r1
  cols c0:c1 …`) to every window so you can confirm at a glance they all match.

`get_roi_bounds()` returns `(row_start, row_end, col_start, col_end)` — the exact
slice used during a scan.

---

## 9. Measurement windows

- **Live View** (`sub_live_lw.py`) — continuous stream for alignment and ROI
  definition. Click a pixel for row/column/time profiles; draw the ROI for a
  time-trace of its mean. Also hosts a delay-stage panel. **This is where the
  shared ROI is set.**
- **ΔT Live + Stage** (`sub_deltat_lw.py`) — live ΔT/T image with manual delay
  control (absolute move + fs step buttons).
- **Pump–Probe Scan** (`sub_pumpprobe_lw.py`) — step the delay stage over up to
  three user-defined fs intervals; plot signal vs delay; save a datacube.
- **TWINS FTIR** (`sub_twins_lw.py`) — scan the GEMINI interferometer at fixed
  delay → interferogram → **power spectrum** (magnitude; no phasing). Static
  hyperspectral imaging.
- **TWINS Pump–Probe** (`sub_twins_pumpprobe_lw.py`) — nested scan: for each pump
  delay, run a GEMINI interferogram → **phase-corrected** spectrum, building a
  hyperspectral map (delay × wavelength). The most processing-heavy mode; see
  [§10](#10-twins-spectral-processing).
- **K-Space Hyperspectral** (`sub_kspace_lw.py`) — like TWINS but keeps **every
  ROI pixel**: builds a datacube `(N_positions, h, w)`, then per-pixel DFT →
  spectral cube `(n_freq, h, w)`. Wavelength slider → spatial map; click a pixel →
  its spectrum.
- **Lock-In Pump–Probe** (`sub_lockin_pumpprobe_lw.py`) — single-wavelength
  time-resolved scan through the SR865A lock-in instead of the camera.

---

## 10. TWINS spectral processing

Implemented in `SpectrumProcessor` (one in `sub_twins_lw.py`, a phase-correcting
one in `sub_twins_pumpprobe_lw.py`). An interferogram is a 1-D trace of the ROI
signal vs GEMINI stage position; the spectrum is its Fourier transform mapped to
wavelength via calibration.

### 10.1 Pipeline (`compute_complex_spectrum`)
1. **Baseline removal** — subtract a moving-average baseline.
2. **Motor-jitter position calibration** — `calibrate_position_axis()` (§10.2).
3. **Centerburst (ZPD) detection** — §10.3.
4. **Optional symmetrization** — mirror the long side about ZPD (for asymmetric
   scans), so a single-sided interferogram becomes double-sided.
5. **Apodization** — Gaussian window centered on ZPD (NIREOS formula).
6. **DFT** — an *explicit* DFT evaluated at frequencies chosen inside the
   `[wl_start, wl_stop]` band and mapped to wavelength by the spectral calibration
   (not a uniform FFT).

### 10.2 Motor-jitter position calibration
The GEMINI's nominal position grid has a small reproducible nonlinearity. Using a
stored reference interferogram (`parameters_int.txt`), `get_real_position_axis()`
recovers the true axis via the analytic-signal phase, and
`calibrate_position_axis()` applies it to every scan. Falls back to the raw axis
if the calibration file is absent.

### 10.3 Centerburst (ZPD) detection
ZPD (zero optical-path difference) is the centerburst — it sets the apodization
center and the linear phase term ($x_0$ in [§10.4](#104-phase-correction-twins-pumpprobe-autophase)).
Detection (`_find_center`) takes the maximum of the **analytic-signal envelope**
$\big|\mathcal H[\tilde I](x)\big|$ (the Hilbert transform magnitude) restricted to
a window around an expected ZPD $x_c$ of half-width $w$:

$$x_{\text{ZPD}} = \underset{|x - x_c|\le w}{\arg\max}\;
   \big|\mathcal H[\tilde I](x)\big|.$$

- The **envelope** is sign-agnostic and robust to which fringe is the local peak,
  and to symmetric/asymmetric scans (the centerburst is the global envelope max
  regardless of where it sits in the array).
- The **window** $[x_c-w,\,x_c+w]$ absorbs the small run-to-run drift while
  rejecting spurious maxima elsewhere in the trace.
- Falls back to a plain `|signal|` argmax if SciPy/positions are unavailable.

Editable in the TWINS Pump–Probe **FFT Settings**: **`ZPD (mm)`** (default
`24.33`) and **`ZPD win (±mm)`** (default `0.1`; set `0` to search the whole
scan). Detected on the reference and reused for all data and all quantities.

### 10.4 Phase correction (TWINS pump–probe "autophase")

**Why phasing is needed.** An interferogram $I(x)$ (signal vs. calibrated stage
position $x$) is real, so its Fourier transform is Hermitian and carries a
spectral phase. After baseline removal and apodization, the windowed
interferogram is

$$\tilde I(x_j) = W(x_j)\,\big[I(x_j) - \bar I(x_j)\big],$$

and the complex spectrum is the (non-uniform) discrete Fourier transform that the
code evaluates directly on the calibrated axis — `compute_complex_spectrum`:

$$S(\nu) = \sum_j \tilde I(x_j)\, e^{-2\pi i\,\nu x_j}\,\Delta x_j
        = |S(\nu)|\,e^{\,i\phi(\nu)},$$

where the spatial frequencies $\nu$ are mapped to wavelength $\lambda$ by the
spectral calibration ([§10.1](#101-pipeline-compute_complex_spectrum)) and
$\Delta x_j$ is the local sample spacing.

For an interferogram that is **perfectly symmetric about ZPD**, $S(\nu)$ is real
($\phi\equiv 0$). Real measurements are not: two effects tilt $\phi(\nu)$ —

1. **Instrumental dispersion** of the birefringent interferometer → a smooth,
   slowly varying $\phi(\nu)$.
2. **A ZPD sampling offset** $x_0$ (the centerburst rarely lands exactly on a
   sample) → a linear phase ramp $\phi(\nu)\approx 2\pi\nu x_0$.

Taking $|S(\nu)|$ (a power spectrum, `compute_spectrum`) throws away the sign and
is fine for static FTIR, but pump–probe needs the **signed** $\Delta T/T$. So we
must rotate each spectral component back onto the real axis.

**Reference-based correction.** The pump-off (odd) interferogram carries the
*same* instrument dispersion and ZPD offset as the data, so it estimates the
phase to remove. From the reference complex spectrum $S_{\text{ref}}(\nu)$
(`compute_phase_correction`):

$$\phi_{\text{ref}}(\nu) = \mathrm{unwrap}\,\arg S_{\text{ref}}(\nu).$$

Because $\arg$ is noisy where $|S_{\text{ref}}|$ is small (band edges, absorption
dips), the raw per-bin phase is **not** used directly. Instead it is replaced by a
magnitude-weighted least-squares polynomial fit of degree $p=$ `PHASE_FIT_ORDER`
($=5$), exploiting that the true instrument phase is smooth in $\nu$:

$$\phi_{\text{corr}} = \arg\min_{\deg P \le p}\;
   \sum_\nu \big|S_{\text{ref}}(\nu)\big|^{2}\,
   \big[\phi_{\text{ref}}(\nu) - P(\nu)\big]^{2}.$$

**Applying it.** Each data spectrum is de-rotated and the in-phase part kept
(`compute_phased_spectrum`):

$$S_{\text{corr}}(\nu) = S_{\text{data}}(\nu)\,e^{-i\,\phi_{\text{corr}}(\nu)},
  \qquad
  \Delta T/T(\lambda) = \mathrm{Re}\big\{S_{\text{corr}}(\nu)\big\}.$$

The discarded $\mathrm{Im}\{S_{\text{corr}}\}$ is the residual dispersive
component. The reference and every data scan share the same ZPD index and the same
frequency grid (`pad_length`), so $\phi_{\text{corr}}$ aligns bin-for-bin with
$S_{\text{data}}$.

**Sign and special cases.** `Invert Polarity` applies an overall
$\Delta T/T \to -\Delta T/T$ if the reference phase is inverted relative to the
transient. A transmission like $T_\text{on}=$ even is a real, non-negative
quantity with no meaningful interferometric phase; with a good reference,
$\mathrm{Re}\{S\,e^{-i\phi_{\text{corr}}}\}\approx |S|$, so passing it
through the same pipeline is harmless (it is **not** specially phased).

### 10.5 Spectral points (`n_points`) — smoothness, not resolution
Zero-filling **interpolates** the spectrum (a smoother curve); it does **not** add
resolution (that is fixed by the scan length and OPD range). The **`Pts`** field:

- **`Auto`** (default) = `ZEROFILL_FACTOR (8) × gemini steps`, clamped to
  `[512, 4096]` — smooth and scaling with the scan.
- A manual value (>0) overrides.

For the per-quantity maps, the displayed quantity's spectrum is reused rather than
recomputed (avoids an identical DFT).

---

## 11. Data and file formats

Saves go under `SaveConfig().base_dir` (default `D:\pumpprobedata`) in dated
folders `YYYY/MM/DD/`. NumPy `.npz` bundles are used so all four quantities are
kept and no re-analysis is needed.

**Image/live saves** (`<sample>_<mode>_<HHMMSS>.npz`) contain:
`Ton`, `Tavg`, `DT`, `DT_T`, plus `raw_odd`, `raw_even`. (If no odd/even pair has
been captured yet, a single-image `.npy` fallback is written.)

**TWINS Pump–Probe** writes, per delay step, `…_step<NNN>_<fs>fs.npz`
(`delay_fs`, `wavelengths`, `spectrum`, `delta_t`, `interferogram`,
`gemini_positions` + calibrated, and the selected per-quantity datacubes), and a
final `…_FINAL.npz`:

- `time_points_fs`, `wavelengths`, `hyperspectral_map` (displayed quantity),
- `reference_wavelengths`, `reference_spectrum`, `zero_mm`,
- **one hyperspectral map per checked quantity**: `hyperspectral_map_Ton`,
  `hyperspectral_map_Tavg`, `hyperspectral_map_DT`, `hyperspectral_map_DT_T`
  (same delay × wavelength axes), so any quantity can be plotted directly:

```python
import numpy as np
d = np.load("…_twins_pp_HHMMSS_FINAL.npz")
ton_map = d["hyperspectral_map_Ton"]   # (time, wavelength), ready to plot
```

---

## 12. Configuration and persisted settings

- **Save directory** — `SaveConfig` (set in the launcher; default `D:\pumpprobedata`).
- **Shared ROI** — `live_roi_state.json` next to the code (auto-managed).
- **TWINS scan/FFT settings** — persisted via `QSettings('Polimi','HybridCamera')`:
  gemini start/stop/steps, wavelength band, apodization, `n_points`, invert/
  symmetrize, **expected ZPD + window**, sample name, settle time.

---

## 13. Troubleshooting

- **LabVIEW "cannot find subVI" after moving the folder** — LabVIEW stores subVI
  links by absolute path. Open `Experiment_manager.vi` *in the LabVIEW editor*,
  point the search dialog to the subVI in this folder, then **File → Save**.
  Alternatively add the project folder to *Tools → Options → Paths → VI Search
  Path*.
- **`pywin32 not installed` / RPC error** — `pip install pywin32`. If hooking a
  new LabVIEW instance fails, open and run `Experiment_manager.vi` manually, then
  Start Manager again.
- **Camera won't initialise** — confirm the `Camera/pt_2dmct1_*.dll` SDK is
  reachable by the VIs and the camera is powered/enumerated.
- **TWINS stage red LED + dialog** — the SmarAct `SCU3DControl.dll` isn't loadable
  or the stage isn't connected (it intentionally does **not** simulate). Fix the
  driver/cabling and retry.
- **Delay stage shows a blue LED** — simulated connection (no Thorlabs hardware
  found); fine for UI testing, not for real moves.
- **Flat/garbage spectra** — acquire a **Reference** before a TWINS pump–probe
  scan (needed for phasing), check the **ZPD** field is near the real burst, and
  verify the calibration files loaded (console `[OK] Loaded … calibration`).
- **Low disk space** — saves write four quantities (≈4× a single image); point the
  save directory at a drive with room.

---

## 14. Developer notes

- Drivers and `ROIState`/`SaveConfig` are **singletons** — one instance per
  process, shared across all windows.
- Missing optional drivers degrade gracefully (`[WARN]` + disabled feature),
  so the camera-only path runs without the stage/lock-in packages.
- The two `SpectrumProcessor` classes are intentionally separate: `sub_twins_lw`
  computes a power spectrum (magnitude); `sub_twins_pumpprobe_lw` adds reference
  phase correction for signed ΔT/T.
- Hyperspectral maps for all selected quantities are computed online through one
  shared pipeline (`_spectrum_from_ifg`), so saved maps need no offline
  re-analysis.
- **Saved data is self-describing**: every `.npz`/`.npy` carries a `meta` JSON
  string (`acq_metadata.meta_json`) with the sample, ROI, scan/FFT/ZPD/phase
  settings, background state, and timestamp. Read it back with
  `json.loads(str(np.load(path)["meta"]))`.
- **Tests** for the numeric core live in `tests/test_processing.py` (metadata,
  ROI bounds, calibration, ZPD detection, phase smoothing, spectrum transform).
  Run `python tests/test_processing.py` (standalone, no deps) or `pytest -q`.
