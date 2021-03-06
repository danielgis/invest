# coding=UTF-8
import sys
import os
import itertools
import glob
from PyInstaller.compat import is_win, is_darwin

# Global Variables
current_dir = os.getcwd()  # assume we're building from the project root
block_cipher = None
exename = 'invest'


kwargs = {
    'hookspath': [os.path.join(current_dir, 'exe', 'hooks')],
    'excludes': None,
    'pathex': sys.path,
    'runtime_hooks': [os.path.join(current_dir, 'exe', 'hooks', 'rthook.py')],
    'hiddenimports': [
        'natcap',
        'natcap.invest',
        'natcap.invest.ui.launcher',
        'yaml',
        'distutils',
        'distutils.dist',
        'rtree',  # mac builds aren't picking up rtree by default.
    ],
    'datas': [('qt.conf', '.')],
    'cipher': block_cipher,
}

cli_file = os.path.join(current_dir, 'src', 'natcap', 'invest', 'cli.py')
a = Analysis([cli_file], **kwargs)

# Compress pyc and pyo Files into ZlibArchive Objects
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Create the executable file.
if is_darwin:
    # Avoid shapely and matplotlib dylib collision with GDAL dylibs.
    excluded_dylibs = set(['libgeos_c.1.dylib', 'libpng16.16.dylib'])
    a.binaries = [x for x in a.binaries if x[0] not in excluded_dylibs]

    # add gdal dynamic libraries from homebrew
    a.binaries += [('geos_c.dll', '/usr/local/lib/libgeos_c.dylib', 'BINARY')]
    a.binaries += [
        (os.path.basename(name), name, 'BINARY') for name in
         itertools.chain(
            glob.glob('/usr/local/lib/libgeos*.dylib'),
            glob.glob('/usr/local/lib/libgeotiff*.dylib'),
            glob.glob('/usr/local/lib/libpng*.dylib'),
            glob.glob('/usr/local/lib/libspatialindex*.dylib')
        )]
elif is_win:
    # Adapted from
    # https://shanetully.com/2013/08/cross-platform-deployment-of-python-applications-with-pyinstaller/
    # Supposed to gather the mscvr/p DLLs from the local system before
    # packaging.  Skirts the issue of us needing to keep them under version
    # control.
    a.binaries += [
        ('msvcp90.dll', 'C:\\Windows\\System32\\msvcp90.dll', 'BINARY'),
        ('msvcr90.dll', 'C:\\Windows\\System32\\msvcr90.dll', 'BINARY')
    ]

    # .exe extension is required if we're on windows.
    exename += '.exe'

exe = EXE(
    pyz,
    a.scripts,
    name=exename,
    exclude_binaries=True,
    debug=False,
    strip=False,
    upx=False,
    console=True)

# Collect Files into Distributable Folder/File
dist = COLLECT(
        exe,
        a.binaries,
        a.zipfiles,
        a.datas,
        name="invest",  # name of the output folder
        strip=False,
        upx=False)
