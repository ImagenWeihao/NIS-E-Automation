# SpheroidPA.spec  --  PyInstaller build spec
# Lives in SpheroidPA\installer\; sources are in the parent directory.
# Run build_installer.bat -- do not run pyinstaller directly.

import os as _os
_SRC = _os.path.abspath(_os.path.join(SPECPATH, '..'))

block_cipher = None

a = Analysis(
    [_os.path.join(_SRC, 'spheroid_pa_gui.py')],
    pathex=[_SRC],
    binaries=[],
    datas=[
        # NIS-E macros -- copied next to the exe so the operator can access them
        (_os.path.join(_SRC, 'nis_macro_auto_capture.mac'), '.'),
        (_os.path.join(_SRC, 'nis_macro_z_autofocus.mac'), '.'),
        (_os.path.join(_SRC, 'nis_macro_20x_capture.mac'), '.'),
    ],
    hiddenimports=[
        # nd2, scipy, skimage imported inside try/except -- PyInstaller misses them
        'nd2',
        'nd2._util',
        'nd2.structures',
        'scipy.ndimage',
        'scipy.ndimage._filters',
        'scipy.ndimage._measurements',
        'scipy.ndimage._morphology',
        'skimage',
        'skimage.measure',
        'skimage.feature',
        'skimage.segmentation',
        'skimage.util',
        'skimage.filters',
        # matplotlib TkAgg backend is loaded dynamically
        'matplotlib',
        'matplotlib.backends.backend_tkagg',
        'matplotlib.backends._backend_tk',
        # local modules
        'spheroid_pipeline',
        'spheroid_screener',
        'cross_zoom_v2',
        'csv_to_nis_bin',
        'nis_bin_to_csv',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'cellpose',
        'watchdog',
        'torch',
        'tensorflow',
        'llvmlite',
        'numba',
        'pyarrow',
        'pandas',
        'h5py',
        'pydantic',
        'pydantic_core',
        'IPython',
        'jupyter',
        'pytest',
        'Pythonwin',
        'win32com',
        'google',
        'grpc',
        'boto',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='SpheroidPA',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='SpheroidPA',
)
