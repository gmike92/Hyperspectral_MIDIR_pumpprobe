"""
Unit tests for the numeric core — the signal math where a silent error would
corrupt data without crashing the GUI:

  * acquisition metadata serialization (acq_metadata)
  * shared ROI bounds (roi_state)
  * stage-axis / spectral calibration (calibration)
  * the TWINS pump-probe SpectrumProcessor: ZPD detection, phase smoothing,
    and the interferogram -> spectrum transform.

Runs either with pytest:

    pip install pytest && pytest -q

or standalone (no pytest needed):

    python tests/test_processing.py
"""

import os
import sys
import math
import json

import numpy as np

# Make the repo root importable whether run via pytest or directly.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# acq_metadata
# ---------------------------------------------------------------------------
def test_meta_json_roundtrip_and_timestamp():
    from acq_metadata import meta_json
    d = json.loads(meta_json(experiment="x", sample="s", roi_bounds=(1, 2, 3, 4)))
    assert d["experiment"] == "x"
    assert d["sample"] == "s"
    assert d["roi_bounds"] == [1, 2, 3, 4]   # tuple -> list through JSON
    assert "saved_at" in d                   # timestamp always added


def test_meta_json_coerces_nonserializable():
    from acq_metadata import meta_json
    d = json.loads(meta_json(obj=object()))  # not JSON-serializable
    assert isinstance(d["obj"], str)         # coerced, save never fails


# ---------------------------------------------------------------------------
# roi_state
# ---------------------------------------------------------------------------
def test_roi_bounds_math():
    from roi_state import ROIState
    s = ROIState()
    s.pos = [100.0, 50.0]    # [col(x), row(y)] top-left
    s.size = [40.0, 25.0]    # [w, h]
    # bounds are (row_start, row_end, col_start, col_end)
    assert s.get_roi_bounds() == (50, 75, 100, 140)


# ---------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------
def test_get_real_position_axis_normalized_and_same_length():
    from calibration import get_real_position_axis
    x = np.linspace(0.0, 1.0, 512)
    ref = np.cos(2 * np.pi * 20 * x)
    axis = get_real_position_axis(ref)
    assert axis.shape == ref.shape
    assert math.isclose(float(axis.min()), 0.0, abs_tol=1e-9)
    assert math.isclose(float(axis.max()), 1.0, abs_tol=1e-9)


def test_calibrate_position_axis_preserves_length():
    from calibration import calibrate_position_axis
    pos = np.linspace(23.8, 24.8, 120)
    out = np.asarray(calibrate_position_axis(pos))
    assert out.shape[0] == pos.shape[0]


# ---------------------------------------------------------------------------
# SpectrumProcessor (TWINS pump-probe)
# ---------------------------------------------------------------------------
def _processor():
    from sub_twins_pumpprobe_lw import SpectrumProcessor
    return SpectrumProcessor()


def _burst(pos, zpd=24.33, fringe=0.03, width=0.05):
    """A realistic oscillatory centerburst at `zpd`."""
    return np.cos(2 * np.pi * (pos - zpd) / fringe) * np.exp(-((pos - zpd) / width) ** 2)


def test_find_center_window_rejects_spike():
    proc = _processor()
    pos = np.linspace(23.8, 24.8, 400)
    sig = _burst(pos)
    sig[50] = 5.0                       # rogue spike at ~23.93 mm, outside window
    proc.zpd_expected, proc.zpd_window = 24.33, 0.1
    idx = proc._find_center(sig, pos)
    assert abs(pos[idx] - 24.33) < 0.02  # lands on the true burst, not the spike


def test_find_center_envelope_without_window():
    proc = _processor()
    pos = np.linspace(23.8, 24.8, 400)
    sig = _burst(pos)
    proc.zpd_expected, proc.zpd_window = None, None
    idx = proc._find_center(sig, pos)
    assert abs(pos[idx] - 24.33) < 0.03


def test_smooth_phase_reduces_noise():
    proc = _processor()
    n = 400
    x = np.linspace(-1.0, 1.0, n)
    true = 0.8 * x + 0.3 * x ** 2 - 0.5 * x ** 3
    mag = np.exp(-(x / 0.6) ** 2) + 0.02            # weak at the edges
    rng = np.random.default_rng(0)
    raw = true + rng.normal(0.0, 0.5, n) / np.sqrt(mag)  # noisy where weak
    fit = proc._smooth_phase(raw, mag)
    rmse_raw = float(np.sqrt(np.mean((raw - true) ** 2)))
    rmse_fit = float(np.sqrt(np.mean((fit - true) ** 2)))
    assert rmse_fit < rmse_raw                      # fit is closer to the truth


def test_compute_spectrum_is_valid_magnitude():
    proc = _processor()
    pos = np.linspace(23.8, 24.8, 400)
    sig = np.cos(2 * np.pi * (pos - 24.3) * 50.0) * np.exp(-((pos - 24.3) / 0.1) ** 2)
    wl, spec = proc.compute_spectrum(pos, sig, n_points=1024, wl_start=8.0, wl_stop=14.0)
    assert wl is not None and spec is not None
    assert len(wl) == len(spec) == 1024
    assert np.all(np.isfinite(spec))
    assert np.all(spec >= 0.0)              # compute_spectrum returns a magnitude
    assert spec.max() > spec.min()         # not a flat/degenerate spectrum
    # wavelengths are returned ascending across the requested band
    assert wl[0] < wl[-1]


# ---------------------------------------------------------------------------
# Standalone runner (no pytest required)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except Exception as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e!r}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
