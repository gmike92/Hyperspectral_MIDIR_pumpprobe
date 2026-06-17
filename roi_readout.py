"""
Small shared UI helper: a compact, auto-refreshing label that shows the current
shared ROI bounds. Dropped into every sub window so the user can verify at a
glance that the ROI is identical across all UIs (it reads the same ROIState).
"""

from pyqtgraph.Qt import QtCore, QtWidgets

from roi_state import ROIState


def add_roi_readout(window, layout):
    """Add a small shared-ROI readout label to `layout`, refreshed every second.

    The label/timer are stored on `window` so they are not garbage-collected.
    Returns the label.
    """
    roi_state = ROIState()
    lbl = QtWidgets.QLabel(roi_state.summary())
    lbl.setWordWrap(True)
    lbl.setStyleSheet("color:#9E9E9E; font-size:10px; font-style:italic; padding:2px;")
    lbl.setToolTip("Shared ROI defined in Live View — identical across all windows.")
    layout.addWidget(lbl)

    timer = QtCore.QTimer(window)
    timer.timeout.connect(lambda: lbl.setText(roi_state.summary()))
    timer.start(1000)

    window._roi_readout_label = lbl
    window._roi_readout_timer = timer
    return lbl
