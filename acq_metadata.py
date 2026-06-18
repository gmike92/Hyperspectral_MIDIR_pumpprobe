"""
Acquisition metadata helper.

Saved data files (.npz) store arrays but not the settings that produced them.
This serializes a small metadata dict to a JSON string so every save is
self-describing (sample, ROI, scan/FFT parameters, background state, timestamp)
and can be re-analyzed later without guessing the acquisition conditions.

Usage:
    from acq_metadata import meta_json
    np.savez(path, DT=..., meta=meta_json(sample=name, roi_bounds=bounds, ...))

Read back:
    import json, numpy as np
    d = np.load(path, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
"""

import json
from datetime import datetime


def meta_json(**fields):
    """Return a JSON string of the given metadata fields.

    A `saved_at` ISO timestamp is always added. Anything not JSON-serializable
    is coerced to str so a save never fails because of metadata.
    """
    base = {"saved_at": datetime.now().isoformat(timespec="seconds")}
    base.update(fields)
    try:
        return json.dumps(base, default=str)
    except Exception:
        # Last-resort: stringify the whole thing so saving still succeeds.
        return json.dumps({k: str(v) for k, v in base.items()})
