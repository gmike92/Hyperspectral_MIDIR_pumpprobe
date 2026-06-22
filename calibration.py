"""
Shared calibration loaders for TWINS / Gemini stage.

Two calibration files live under ./Twins/ASRC calibration/:

    parameters_cal.txt   Spectral calibration: wavelength (µm) <-> reciprocal (1/mm-stage).
                         Used by every spectrum builder to map stage-Fourier
                         frequencies onto real wavelengths.

    parameters_int.txt   Motor-nonlinearity reference: an interferogram of a
                         known-wavelength source vs. nominal stage positions.
                         Lets us recover the *real* position axis (the analytic-
                         signal phase trick from NIREOS get_real_position_axis)
                         and remove the motor's reproducible nonlinearity from
                         every scan.

Both are cached at module load: each file is opened at most once per process.
Every TWINS-style processor (sub_twins_lw, sub_twins_pumpprobe_lw,
sub_kspace_lw) calls into here instead of reloading independently.
"""

from pathlib import Path
import numpy as np


_SPECTRAL_PATHS = [
    Path(r".\Twins\ASRC calibration\parameters_cal.txt"),
    Path(r"C:\Users\mguizzardi\Desktop\Camera python\TWINS FILE\Twins\ASRC calibration\parameters_cal.txt"),
]

_POSITION_PATHS = [
    Path(r".\Twins\ASRC calibration\parameters_int.txt"),
    Path(r"C:\Users\mguizzardi\Desktop\Camera python\TWINS FILE\Twins\ASRC calibration\parameters_int.txt"),
]


# Module-level caches: None = not yet attempted, (a, b) = loaded, (None, None) = tried and failed
_spectral_cache = None
_position_cache = None


def _read_two_row_file(path):
    """Read a tab-separated file with two rows: row 0 = x-axis, row 1 = y-axis."""
    import pandas as pd
    ref = pd.read_csv(path, sep="\t", header=None)
    return (ref.iloc[0].to_numpy(dtype='float64'),
            ref.iloc[1].to_numpy(dtype='float64'))


def get_spectral_calibration():
    """Return (wavelength_cal_um, reciprocal_cal_per_mm) or (None, None)."""
    global _spectral_cache
    if _spectral_cache is not None:
        return _spectral_cache
    for p in _SPECTRAL_PATHS:
        if p.exists():
            try:
                wl, rk = _read_two_row_file(p)
                print(f"[OK] Loaded spectral calibration: {p.name}  "
                      f"({wl.min():.2f}–{wl.max():.2f} µm)")
                _spectral_cache = (wl, rk)
                return _spectral_cache
            except Exception as e:
                print(f"[WARN] Could not read spectral calibration {p}: {e}")
    print("[WARN] parameters_cal.txt not found — using simple 1/freq conversion.")
    _spectral_cache = (None, None)
    return _spectral_cache


def get_position_calibration():
    """Return (position_ref_mm, amplitude_ref) or (None, None).

    These are the nominal-stage positions and the reference-source interferogram
    used by calibrate_position_axis().
    """
    global _position_cache
    if _position_cache is not None:
        return _position_cache
    for p in _POSITION_PATHS:
        if p.exists():
            try:
                pos, amp = _read_two_row_file(p)
                print(f"[OK] Loaded position calibration: {p.name}")
                _position_cache = (pos, amp)
                return _position_cache
            except Exception as e:
                print(f"[WARN] Could not read position calibration {p}: {e}")
    print("[WARN] parameters_int.txt not found — motor jitter correction disabled.")
    _position_cache = (None, None)
    return _position_cache


def calibration_status():
    """Whether the calibrations are actually loaded (vs. falling back to an
    uncalibrated axis). For embedding in saved metadata so a file records that
    its position axis was motor-jitter corrected.

        position_axis_calibrated : motor-nonlinearity (parameters_int.txt) applied
        spectral_calibrated      : wavelength mapping (parameters_cal.txt) applied
    """
    wl, _ = get_spectral_calibration()
    pos, _ = get_position_calibration()
    return {
        "position_axis_calibrated": pos is not None,
        "spectral_calibrated": wl is not None,
    }


def get_real_position_axis(reference):
    """Recover the real position axis from a reference interferogram.

    Analytic-signal phase trick (NIREOS get_real_position_axis):
        FFT → zero negative freqs → iFFT → unwrap(-angle) → normalize to [0,1].
    """
    ref = np.asarray(reference).squeeze()
    fft_ref = np.fft.fft(ref)
    half = int(np.floor(len(ref) / 2) - 1)
    if half > 0:
        fft_ref[:half] = 0.0
    phase = np.unwrap(-np.angle(np.fft.ifft(fft_ref)))
    a, b = float(phase.min()), float(phase.max())
    if b - a == 0:
        return np.linspace(0.0, 1.0, len(ref))
    return (phase - a) / (b - a)


def calibrate_position_axis(position_axis):
    """Apply motor-nonlinearity correction to a nominal stage axis.

    Interpolates the stored reference IFG onto the user's (oversampled) axis,
    runs get_real_position_axis to extract the corrected axis, downsamples
    back to the input length. Falls back to the input unchanged if the
    position calibration isn't loaded.
    """
    pos_ref, amp_ref = get_position_calibration()
    if pos_ref is None or amp_ref is None:
        return np.asarray(position_axis)

    try:
        from scipy import interpolate as _interp
        position_axis = np.asarray(position_axis).squeeze()
        if position_axis.size < 4:
            return position_axis

        d_user = float(np.mean(np.diff(position_axis)))
        d_ref = float(np.mean(np.diff(pos_ref)))
        if d_ref == 0:
            return position_axis
        factor = int(max(1, np.ceil(abs(d_user / d_ref))))

        oversmp_pos = np.linspace(position_axis[0], position_axis[-1],
                                  factor * position_axis.size)
        f_ref = _interp.interp1d(pos_ref, amp_ref, kind='cubic',
                                 bounds_error=False,
                                 fill_value=(amp_ref[0], amp_ref[-1]))
        ref_interp = f_ref(oversmp_pos)

        calib_norm = get_real_position_axis(ref_interp)
        calib_overs = calib_norm * (position_axis[-1] - position_axis[0]) + position_axis[0]
        calibrated = calib_overs[0:-1:factor]
        if calibrated.size >= position_axis.size:
            return calibrated[:position_axis.size]
        pad = np.full(position_axis.size - calibrated.size,
                      calibrated[-1] if calibrated.size else 0.0)
        return np.concatenate([calibrated, pad])
    except Exception as e:
        print(f"[WARN] Position calibration failed, using raw axis: {e}")
        return np.asarray(position_axis)
