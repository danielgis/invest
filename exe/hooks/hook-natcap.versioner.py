from sys import platform as _platform
if _platform == "linux" or _platform == "linux2":
    # linux
    from PyInstaller.hooks.hookutils import collect_data_files, collect_submodules
    datas = collect_data_files('natcap.versioner')
    hiddenimports = collect_submodules('natcap.versioner')
else:
    from PyInstaller.utils.hooks import collect_data_files, collect_submodules
    datas = collect_data_files('natcap.versioner')
    hiddenimports = collect_submodules('natcap.versioner')
