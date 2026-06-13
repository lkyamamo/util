"""
ParaView Python script for shock voxel visualization.

Run inside ParaView's Python shell (Tools → Python Shell) or with pvpython:
    pvpython paraview_script.py trajectory.xdmf [quantity] [output.mp4]

Requires generate_xdmf.py to have been run first to produce the .xdmf file.
"""

import sys

try:
    from paraview.simple import (
        XDMFReader, Clip, Show, Hide, GetActiveViewOrCreate,
        GetAnimationScene, ColorBy, GetColorTransferFunction,
        GetOpacityTransferFunction, RenderAllViews, SaveAnimation,
        ResetCamera,
    )
    import paraview.simple as pv
except ImportError:
    print("paraview.simple not found — run this script with pvpython or inside ParaView.")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
XDMF_FILE = sys.argv[1] if len(sys.argv) > 1 else "trajectory.xdmf"
QUANTITY   = sys.argv[2] if len(sys.argv) > 2 else "density"
OUTPUT     = sys.argv[3] if len(sys.argv) > 3 else None

CMAPS = {
    'density':         'viridis',
    'pressure':        'Cool to Warm',
    'virial_pressure': 'Cool to Warm',
    'temperature':     'Plasma (matplotlib)',
    'avg_speed':       'viridis',
    'avg_O_speed':     'viridis',
    'voxel_type':      'Set1',
}


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
reader = XDMFReader(FileNames=[XDMF_FILE])
reader.CellArrayStatus = [
    'density', 'pressure', 'virial_pressure',
    'temperature', 'avg_speed', 'avg_O_speed', 'voxel_type',
]


# ---------------------------------------------------------------------------
# View setup
# ---------------------------------------------------------------------------
view = GetActiveViewOrCreate('RenderView')
view.Background = [0.1, 0.1, 0.18]    # dark background matching PyVista version


# ---------------------------------------------------------------------------
# Volume representation
# ---------------------------------------------------------------------------
display = Show(reader, view)
display.Representation = 'Volume'
ColorBy(display, ('CELLS', QUANTITY))
display.RescaleTransferFunctionToDataRange(False, True)

cmap = CMAPS.get(QUANTITY, 'viridis')
ctf  = GetColorTransferFunction(QUANTITY)
ctf.ApplyPreset(cmap, True)

# Linear opacity (mirrors opacity='linear' in PyVista)
otf = GetOpacityTransferFunction(QUANTITY)
data_range = display.GetDataInformation().GetCellDataInformation() \
                     .GetArrayInformation(QUANTITY).GetComponentRange(0)
otf.Points = [data_range[0], 0.0, 0.5, 0.0,
              data_range[1], 1.0, 0.5, 0.0]

ResetCamera(view)
RenderAllViews()


# ---------------------------------------------------------------------------
# Optional: clip plane  (uncomment and set values to activate)
# ---------------------------------------------------------------------------
# clip = Clip(Input=reader)
# clip.ClipType = 'Plane'
# clip.ClipType.Normal = [1.0, 1.0, 1.0]
# clip.ClipType.Origin = [25.0, 25.0, 25.0]
# clip.InsideOut = False
# Hide(reader, view)
# display = Show(clip, view)
# display.Representation = 'Volume'
# ColorBy(display, ('CELLS', QUANTITY))
# display.RescaleTransferFunctionToDataRange()
# RenderAllViews()


# ---------------------------------------------------------------------------
# Animation / video export
# ---------------------------------------------------------------------------
scene = GetAnimationScene()
scene.PlayMode      = 'Sequence'
scene.StartTime     = 0
scene.EndTime       = reader.TimestepValues[-1] if reader.TimestepValues else 0
scene.NumberOfFrames = len(reader.TimestepValues) if reader.TimestepValues else 1

if OUTPUT:
    print(f"Saving animation → {OUTPUT}")
    SaveAnimation(
        OUTPUT, view,
        FrameRate=10,
        ImageResolution=[1920, 1080],
        OverrideColorPalette='',
    )
    print("Done.")
else:
    print("No output path given — interactive mode.")
    print("Use File → Save Animation in the GUI to export video.")
