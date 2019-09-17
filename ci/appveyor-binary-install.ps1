# This is an install script needed for our binary build on AppVeyor.
#
# It assumes that the applications choco, wget, 7zip, unzip and curl are available
# and on the PATH.  It also requires that the environment variable
#

choco install make wget vcredist140 pandoc zip

# The binary build requires the shapely DLL to be named something specific.
# /B copies the file as a binary file.
cmd.exe --% /c copy /B %PYTHON%\Lib\site-packages\shapely\DLLs\geos_c.dll %PYTHON%\Lib\site-packages\shapely\DLLs\geos.dll

# Download and install NSIS plugins to their correct places.
cmd.exe --% /c wget -nv https://storage.googleapis.com/natcap-build-dependencies/windows/Inetc.zip
cmd.exe --% /c wget -nv https://storage.googleapis.com/natcap-build-dependencies/windows/Nsisunz.zip
cmd.exe --% /c wget -nv https://storage.googleapis.com/natcap-build-dependencies/windows/NsProcess.zip
cmd.exe --% /c 7z e NsProcess.zip -o"C:\Program Files (x86)\NSIS\Plugins\x86-ansi" Plugin\nsProcess.dll
cmd.exe --% /c 7z e NsProcess.zip -o"C:\Program Files (x86)\NSIS\Include" Include\nsProcess.nsh
cmd.exe --% /c 7z e Inetc.zip -o"C:\Program Files (x86)\NSIS\Plugins\x86-ansi" Plugins\x86-ansi\INetC.dll
cmd.exe --% /c 7z e Nsisunz.zip -o"C:\Program Files (x86)\NSIS\Plugins\x86-ansi" nsisunz\Release\nsisunz.dll
