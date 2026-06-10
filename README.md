# LabVIEW Hybrid Camera System

A hybrid Python/LabVIEW control software suite designed for time-resolved optical spectroscopy, featuring modular architecture for various types of data acquisition. The software integrates high-level Python Graphical User Interfaces (built with PyQtGraph/PyQt6) with underlying LabVIEW VIs (Virtual Instruments), Thorlabs delay stages, and Stanford Research Lock-In amplifiers.

## Overview and Architecture

The system uses a "puppeteer" design pattern where Python handles the UI, scan logic, data plotting, and file I/O, while delegating low-level camera interactions and initial signal processing to an actively running LabVIEW VI (`Experiment_manager.vi`).

### Core Components
- **`main_launcher_lw.py`**: The primary entry point for launching the software. It provides a global control panel to start the LabVIEW Manager, connect to hardware stages, and open various sub-module windows.
- **`labview_manager.py` & `hybrid_camera_server.py`**: These scripts contain `LabVIEWManager`, which uses ActiveX (`win32com.client`) to command the running `Experiment_manager.vi` in the background. It reads data (Transmission $T$ and Differential Transmission $\Delta T/T$) computed continuously by LabVIEW.
- **`stage_delay.py`**: Python controller for the optical delay stage (often Thorlabs hardware). Supports moving, homing, and simulated execution.
- **`stage_driver.py`**: Python controller for the "Twins" interferometer stage.
- **`driver_lockin.py`**: Driver for the Stanford Research SR865A Lock-In amplifier via QCoDeS, enabling non-camera measurements.

## Sub-Modules (Measurement Modes)

The system is separated into distinct sub-windows, each specialized for a specific type of experiment:

1. **Live View (`sub_live_lw.py`)** 📹
   Provides a direct, real-time feed of the camera sensor for alignment, debugging, and general beam monitoring.
   
2. **Pump-Probe Scan (`sub_pumpprobe_lw.py`)** 🔬
   Conducts standard Pump-Probe measurements using an optical delay line and the camera detector. Features multiple user-defined temporal intervals tailored to track ultra-fast dynamics at various resolutions.

3. **$\Delta T/T$ Live + Stage (`sub_deltat_lw.py`)**
   A dedicated module for visualizing dynamic Differential Transmission ($\Delta T/T$) signals in real time while manually or automatically stepping the delay stage.

4. **Twins FTIR (`sub_twins_lw.py`)** 🌊
   Controls the Twins optical interferometer stage to acquire static interferograms and perform Fourier Transform Infrared (FTIR) spectroscopy measurements on the sample without the pump.

5. **Twins Pump-Probe (`sub_twins_pumpprobe_lw.py`)** 🌈
   A multi-dimensional spectroscopy mode matching the Pump-Probe delay stage with the Twins short-delay interferometer, designed for tasks like two-dimensional electronic spectroscopy (2DES).

6. **K-Space Hyperspectral (`sub_kspace_lw.py`)**
   Dedicated module focused on acquiring spatially-resolved interference signals across the camera's active pixels, reconstructing K-space datasets.

7. **Lock-In Pump-Probe (`sub_lockin_pumpprobe_lw.py`)** 🔬
   Performs Pump-Probe measurements acquiring analog signals through a Lock-In amplifier rather than utilizing the camera. Employs `driver_lockin.py` to communicate synchronously with stage movement.

## Installation and Requirements

1. **Python version**: Python 3.8+ (64-bit recommended)
2. **LabVIEW Environment**: A compatible LabVIEW installation containing `Experiment_manager.vi`.
3. **Python Dependencies**:
   Install via `pip`:
   ```bash
   pip install numpy h5py pyqtgraph pyqt6 pywin32 qcodes
   ```
   *(Note: `pywin32` is critical for Python-LabVIEW COM communications).*

## Usage Instructions

1. **Prepare LabVIEW**: Ensure `Experiment_manager.vi` is accessible and properly configured. If LabVIEW is already open, the script will hook into the active instance.
2. **Start the Launcher**:
   ```bash
   python main_launcher_lw.py
   ```
3. **Connect Hardware**:
   - In the launcher, click **"Start Manager"** to initialize communication with LabVIEW.
   - Click **"Init Camera"** to open the sensor inside the VI.
   - Connect the required functional stages (Delay stage and/or Twins stage) under "Stage Controls".
4. **Launch an Experiment**:
   - Once the necessary LEDs turn green (indicating positive connection status), launch the desired operational window from the "Sub-Windows" panel.
   - Each module handles its own data streaming, graph updates, and HDF5 / CSV data saving.

## Troubleshooting

- **LabVIEW Connection Fails ("pywin32 not installed" or RPC error)**: Make sure `pywin32` is correctly installed. If connecting to a new LabVIEW instance fails, manually open `Experiment_manager.vi`, run it, and retry connecting Python.
- **Instrument Driver Missing (`stage_delay.py`, `driver_lockin.py`)**: If physical hardware is missing, the code gracefully falls back or attempts simulation modes depending on the implementation.
- **Lock-In Timeouts**: If the Lock-In amplifier hangs, `driver_lockin.py` forces a Visa handle reset and instrument trace deletion. Verify the USB resource string in `list_visa.py`.
