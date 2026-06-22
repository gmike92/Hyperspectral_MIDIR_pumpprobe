"""
Shared ROI / pixel selection — the single source of truth for the region the
user draws in Live View and that every scan window (k-space, pump-probe, TWINS,
TWINS pump-probe) reads back.

The ROI is defined once in Live View and used everywhere else. A reset or
desynced ROI changes the measured intensity, so it is kept in one place (this
process-wide singleton) and persisted to disk so it survives Live View
close/reopen and full app restarts.

Live View *writes* this store (via update_from) whenever the user moves/resizes
the ROI, toggles ROI/Pixel mode, or clicks a pixel. The scan windows *read* it
(get_roi_bounds / sel_row / sel_col), independent of whether Live View is open.
"""

import os
import json

# Stored next to the modules so every window resolves the same file.
ROI_STATE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "live_roi_state.json")


class ROIState:
    """Process-wide singleton holding the current ROI rectangle, the ROI/Pixel
    toggle, and the selected pixel."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._init_state()
        return cls._instance

    def _init_state(self):
        # Defaults match the historical Live View ROI ([x, y] / [w, h]).
        self.pos = [20.0, 20.0]    # top-left [col (x), row (y)]
        self.size = [30.0, 30.0]   # [width, height]
        self.use_roi = True        # True = ROI mode, False = single-pixel mode
        self.sel_row = None
        self.sel_col = None
        self.load()

    # ------------------------------------------------------------------ #
    # Persistence
    # ------------------------------------------------------------------ #
    def load(self):
        """Reload from disk (call once at construction)."""
        try:
            if os.path.exists(ROI_STATE_FILE):
                with open(ROI_STATE_FILE, "r") as f:
                    d = json.load(f)
                self.pos = d.get("roi_pos", self.pos)
                self.size = d.get("roi_size", self.size)
                self.use_roi = d.get("use_roi", self.use_roi)
                self.sel_row = d.get("sel_row", self.sel_row)
                self.sel_col = d.get("sel_col", self.sel_col)
        except Exception as e:
            print(f"[ROIState] load failed: {e}")

    def save(self):
        """Persist to disk so the selection survives close/reopen and restart."""
        try:
            with open(ROI_STATE_FILE, "w") as f:
                json.dump({
                    "roi_pos": [float(self.pos[0]), float(self.pos[1])],
                    "roi_size": [float(self.size[0]), float(self.size[1])],
                    "use_roi": bool(self.use_roi),
                    "sel_row": self.sel_row,
                    "sel_col": self.sel_col,
                }, f)
        except Exception as e:
            print(f"[ROIState] save failed: {e}")

    # ------------------------------------------------------------------ #
    # Write API (used by Live View)
    # ------------------------------------------------------------------ #
    def update_from(self, pos, size, use_roi, sel_row, sel_col):
        """Update the canonical ROI from the Live View widget and persist it."""
        self.pos = [float(pos[0]), float(pos[1])]
        self.size = [float(size[0]), float(size[1])]
        self.use_roi = bool(use_roi)
        self.sel_row = sel_row
        self.sel_col = sel_col
        self.save()

    # ------------------------------------------------------------------ #
    # Read API (used by every scan window — mirrors the old LiveViewWindow API)
    # ------------------------------------------------------------------ #
    def get_roi_bounds(self):
        """Return current ROI bounds as (row_start, row_end, col_start, col_end)."""
        x0 = max(0, int(self.pos[0]))
        y0 = max(0, int(self.pos[1]))
        rw = max(1, int(self.size[0]))
        rh = max(1, int(self.size[1]))
        return (y0, y0 + rh, x0, x0 + rw)

    def get_roi_shape(self):
        """ROI size in pixels as (n_rows, n_cols) = (height, width)."""
        r0, r1, c0, c1 = self.get_roi_bounds()
        return (r1 - r0, c1 - c0)

    def summary(self):
        """Compact one-line description, shown in every window so the user can
        verify at a glance that the ROI is identical across all sub UIs."""
        r0, r1, c0, c1 = self.get_roi_bounds()
        if self.use_roi:
            return f"Shared ROI  rows {r0}:{r1}  cols {c0}:{c1}  ({r1 - r0}×{c1 - c0} px)"
        return f"Shared ROI  pixel ({self.sel_row}, {self.sel_col})  [Pixel mode]"
