; Variables needed at the command line:
; VERSION         - the version of InVEST we're building (example: 3.4.5)
; VERSION_DISK    - the windows-safe version of InVEST we're building
;                   This needs to be a valid filename on windows (no : , etc),
;                   but could represent a development build.
; INVEST_3_FOLDER - The local folder of binaries to include.
; SHORT_VERSION   - The short version name.  Usually a tagname such as 'tip',
;                   'default', or 3.4.5.
; ARCHITECTURE    - The architecture we're building for.  Generally this is x86.
; FORKNAME        - The username of the InVEST fork we're building off of.
; DATA_LOCATION   - Where (relative to datportal) the data should be downloaded
;                   from.
;
; NOTE ON INSTALLING SAMPLE DATA:
; ===============================
; There are three ways to install sample data with this installer:
;
; 1) Through the Installer's GUI.
;    This approach requires users to interact with the GUI of the installer,
;    where the user will select the data zipfile he/she would like to have
;    installed as part of the installation. If the user does not have an active
;    internet connection (or if there are problems with a download), an error
;    dialog will be presented for each failed download.
;
; 2) Through the 'Advanced' input on the front pane of the installer.
;    This approach is particularly convenient for users wishing to distribute
;    sample data as a single zipfile with the installer, as might be the case
;    for sysadmins installing on many computers, or Natcappers installing on
;    user's computers at a training.  To make this work, a specially formatted
;    zipfile must be used.  This zipfile may be created with paver by calling:
;
;        $ paver build_data --single-archive
;
;    Alternately, this zipfile may be assembled by hand, so long as the
;    zipfile has all sample data folders at the top level.  Whatever is in the
;    archive will be unzipped to the install directory.
;
;    It's also worth noting that this 'Advanced' install may be used at the
;    command-line, optionally as part of a silent install.  If we assume that
;    the InVEST 3.3.1 installer and the advanced sampledata zipfile are copied
;    to the same directory, and we open a cmd prompt within that same
;    directory:
;
;        > .\InVEST_3.3.1_Setup_x86.exe /S /DATAZIP=%CD%\sampledata.zip
;
;    This will execute the installer silently, and extract the contents of
;    sampledata.zip to the installation directory.
;
; 3) By having the installer and sample data archives in the right places
;    This approach is an alternative to the silent install with the 'advanced'
;    input functionality, and is useful when the user has control over the
;    location of the installer and the sampledata zipfiles on the local
;    computer. The gist is that if the installer finds the sample data zipfile
;    it's looking for in the right place, it'll use that instead of going to
;    the network.
;
;    To use this, the following folder structure must exist:
;
;    some directory/
;        InVEST_<version>_Setup.exe
;        sample_data/
;           Marine.zip
;           Pollination.zip
;           Base_Data.zip
;           <other zipfiles, as desired, downloaded from our website>
!include nsProcess.nsh
!include LogicLib.nsh
; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "InVEST"
!define PRODUCT_VERSION "${VERSION} ${ARCHITECTURE}"
!define PDF_NAME "InVEST_${SHORT_VERSION}_Documentation.pdf"
!define PRODUCT_PUBLISHER "The Natural Capital Project"
!define PRODUCT_WEB_SITE "http://www.naturalcapitalproject.org"
!define MUI_COMPONENTSPAGE_NODESC
!define PACKAGE_NAME "${PRODUCT_NAME} ${PRODUCT_VERSION}"

SetCompressor zlib

; MUI has some graphical files that I want to define, which must be defined
; here before the macros are declared.
;
; NOTES ABOUT GRAPHICS:
; ---------------------
; NSIS is surprisingly picky about the sorts of graphics that can be displayed.
; Here's what I know about these images after a fair amount of
; trial and error:
;  * Image format must be Windows Bitmap (.bmp).
;       * I've used 24-bit ad 32-bit encodings without issue.
;       * 24-bit encodings should be sufficient, and yield ~30% filesize reduction.
;       * If using GIMP, be sure to check the compatibility option marked
;         "Do not write color space information".
;  * Vertical images must have dimensions 164Wx314H.
;       * Within this, the InVEST logo currently has dimensions 130Wx109H.
;  * Horizontal (top) banner must have dimensions 150Wx57H.
;       * Within this, the InVEST logo currently has dimensions 48Wx40H.
;
; GIMP notes: I've had good results with just opening the existing BMPs from
; the repo, inserting a new layer with the InVEST logo, scaling the layer,
; repositioning the logo to perfectly cover the old logo, flattening the
; layers and then exporting as a 24-bit windows bitmap.
!define MUI_WELCOMEFINISHPAGE_BITMAP "InVEST-vertical.bmp"
!define MUI_UNWELCOMEFINISHPAGE_BITMAP "InVEST-vertical.bmp"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP "InVEST-header-wcvi-rocks.bmp"
!define MUI_UNHEADERIMAGE_BITMAP "InVEST-header-wcvi-rocks.bmp"
!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\orange-uninstall.ico"

; MUI 1.67 compatible ------
!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "x64.nsh"
!include "FileFunc.nsh"
!include "nsDialogs.nsh"
!include "WinVer.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "InVEST-2.ico"

; Add an advanced options control for the welcome page.
!define MUI_PAGE_CUSTOMFUNCTION_SHOW AddAdvancedOptions
!define MUI_PAGE_CUSTOMFUNCTION_LEAVE ValidateAdvZipFile
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\..\LICENSE.txt"

!define MUI_PAGE_CUSTOMFUNCTION_PRE SkipComponents
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
;!define MUI_FINISHPAGE_SHOWREADME ${PDF_NAME}
!insertmacro MUI_PAGE_FINISH

; MUI Uninstaller settings---------------
!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

; Language files
!insertmacro MUI_LANGUAGE "English"

; MUI end ------

!define INSTALLER_NAME "InVEST_${FORKNAME}${VERSION_DISK}_${ARCHITECTURE}_Setup.exe"
Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile ${INSTALLER_NAME}
InstallDir "C:\InVEST_${VERSION_DISK}_${ARCHITECTURE}"
ShowInstDetails show
RequestExecutionLevel admin

; This function allows us to test to see if a process is currently running.
; If the process name passed in is actually found, a message box is presented
; and the uninstaller should quit.
!macro CheckProgramRunning process_name
    ${nsProcess::FindProcess} "${process_name}.exe" $R0
    Pop $R0

    StrCmp $R0 603 +3
        MessageBox MB_OK|MB_ICONEXCLAMATION "InVEST is still running.  Please close all InVEST models and try again."
        Abort
!macroend

var AdvCheckbox
var AdvFileField
var AdvZipFile
var LocalDataZipFile
Function AddAdvancedOptions
    ${NSD_CreateCheckBox} 120u -18u 15% 12u "Advanced"
    pop $AdvCheckbox
    ${NSD_OnClick} $AdvCheckbox EnableAdvFileSelect

    ${NSD_CreateFileRequest} 175u -18u 36% 12u $LocalDataZipFile
    pop $AdvFileField
    ShowWindow $AdvFileField 0

    ${NSD_CreateBrowseButton} 300u -18u 5% 12u "..."
    pop $AdvZipFile
    ${NSD_OnClick} $AdvZipFile GetZipFile
    ShowWindow $AdvZipFile 0

    ; if $LocalDataZipFile has a value, check the 'advanced' checkbox by default.
    ${If} $LocalDataZipFile != ""
        ${NSD_Check} $AdvCheckbox
        Call EnableAdvFileSelect
    ${EndIf}
FunctionEnd

Function EnableAdvFileSelect
    ${NSD_GetState} $AdvCheckbox $0
    ShowWindow $AdvFileField $0
    ShowWindow $AdvZipFile $0
FunctionEnd

Function GetZipFile
    nsDialogs::SelectFileDialog "open" "" "Zipfiles *.zip"
    pop $0
    ${GetFileExt} $0 $1
    ${If} $1 != "zip"
        MessageBox MB_OK "File must be a zipfile"
        Abort
    ${EndIf}
    ${NSD_SetText} $AdvFileField $0
    strcpy $LocalDataZipFile $0
FunctionEnd

Function SkipComponents
    ${If} $LocalDataZipFile != ""
        Abort
    ${EndIf}
FunctionEnd

Function ValidateAdvZipFile
    ${NSD_GetText} $AdvFileField $0
    ${If} $0 != ""
        ${GetFileExt} $0 $1
        ${If} $1 != "zip"
            MessageBox MB_OK "File must be a zipfile $1"
            Abort
        ${EndIf}
        IfFileExists $0 +3 0
        MessageBox MB_OK "File not found or not accessible: $0"
        Abort
    ${Else}
        ; Save the value in the advanced filefield as $LocalDataZipFile
        strcpy $LocalDataZipFile $0
    ${EndIf}
FunctionEnd

!define LVM_GETITEMCOUNT 0x1004
!define LVM_GETITEMTEXT 0x102D

Function DumpLog
    Exch $5
    Push $0
    Push $1
    Push $2
    Push $3
    Push $4
    Push $6

    FindWindow $0 "#32770" "" $HWNDPARENT
    GetDlgItem $0 $0 1016
    StrCmp $0 0 exit
    FileOpen $5 $5 "w"
    StrCmp $5 "" exit
        SendMessage $0 ${LVM_GETITEMCOUNT} 0 0 $6
        System::Alloc ${NSIS_MAX_STRLEN}
        Pop $3
        StrCpy $2 0
        System::Call "*(i, i, i, i, i, i, i, i, i) i \
            (0, 0, 0, 0, 0, r3, ${NSIS_MAX_STRLEN}) .r1"
        loop: StrCmp $2 $6 done
            System::Call "User32::SendMessageA(i, i, i, i) i \
            ($0, ${LVM_GETITEMTEXT}, $2, r1)"
            System::Call "*$3(&t${NSIS_MAX_STRLEN} .r4)"
            FileWrite $5 "$4$\r$\n"
            IntOp $2 $2 + 1
            Goto loop
        done:
            FileClose $5
            System::Free $1
            System::Free $3
    exit:
        Pop $6
        Pop $4
        Pop $3
        Pop $2
        Pop $1
        Pop $0
        Exch $5
FunctionEnd

Function Un.onInit
    !insertmacro CheckProgramRunning "invest"
FunctionEnd

Section "InVEST Tools and ArcGIS toolbox" Section_InVEST_Tools
  SetShellVarContext all
  SectionIn RO ;require this section

  !define SMPATH "$SMPROGRAMS\${PACKAGE_NAME}"
  !define INVEST_ICON "$INSTDIR\${INVEST_3_FOLDER}\InVEST-2.ico"
  !define INVEST_DATA "$INSTDIR\${INVEST_3_FOLDER}"
  !define OVERLAP "${SMPATH}\Overlap Analysis"
  !define HRA "${SMPATH}\Habitat Risk Assessment"
  !define COASTALBLUECARBON "${SMPATH}\Coastal Blue Carbon"
  !define FISHERIES "${SMPATH}\Fisheries"
  !define HYDROPOWER "${SMPATH}\Hydropower"

  ; Write the uninstaller to disk
  SetOutPath "$INSTDIR"
  !define UNINSTALL_PATH "$INSTDIR\Uninstall_${VERSION_DISK}.exe"
  writeUninstaller "${UNINSTALL_PATH}"

  ; Create start  menu shortcuts.
  ; These shortcut paths are set in the appropriate places based on the SetShellVarConext flag.
  ; This flag is automatically set based on the MULTIUSER installation mode selected by the user.
  SetOutPath "$INSTDIR\${INVEST_3_FOLDER}"

  CreateDirectory "${SMPATH}"
  CreateShortCut "${SMPATH}\Crop Production (Percentile) (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_crop_production_percentile.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Crop Production (Regression) (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_crop_production_regression.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Scenic Quality (unstable) (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_scenic_quality.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Habitat Quality (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_habitat_quality.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Carbon (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_carbon.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Forest Carbon Edge Effect (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_forest_carbon_edge_effect.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\GLOBIO (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_globio.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Pollination (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_pollination.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Finfish Aquaculture (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_finfish_aquaculture.bat" "" "${INVEST_ICON}"
  CreateDirectory "${OVERLAP}"
  CreateShortCut "${OVERLAP}\Overlap Analysis (Management Zones) (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_overlap_analysis_mz.bat" "" "${INVEST_ICON}"
  CreateShortCut "${OVERLAP}\Overlap Analysis (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_overlap_analysis.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Wave Energy (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_wave_energy.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Wind Energy (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_wind_energy.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Coastal Vulnerability (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_coastal_vulnerability.bat" "" "${INVEST_ICON}"

  CreateDirectory "${COASTALBLUECARBON}"
  CreateShortCut "${COASTALBLUECARBON}\(1) Coastal Blue Carbon Preprocessor (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_coastal_blue_carbon_preprocessor.bat" "" "${INVEST_ICON}"
  CreateShortCut "${COASTALBLUECARBON}\(2) Coastal Blue Carbon (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_coastal_blue_carbon.bat" "" "${INVEST_ICON}"

  CreateDirectory "${FISHERIES}"
  CreateShortCut "${FISHERIES}\(1) Fisheries (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_fisheries.bat" "" "${INVEST_ICON}"
  CreateShortCut "${FISHERIES}\(2) Fisheries Habitat Scenario Tool (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_fisheries_hst.bat" "" "${INVEST_ICON}"

  CreateDirectory "${HRA}"
  CreateShortCut "${HRA}\(1) Habitat Risk Assessment Preprocessor (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_hra_preprocessor.bat" "" "${INVEST_ICON}"
  CreateShortCut "${HRA}\(2) Habitat Risk Assessment (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_hra.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\SDR (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_sdr.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\NDR (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_ndr.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Scenario Generator: Rule Based (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_scenario_generator.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Scenario Generator: Proximity Based (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_scenario_gen_proximity.bat" "" "${INVEST_ICON}"

  CreateShortCut "${SMPATH}\Water Yield (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_hydropower_water_yield.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\Seasonal Water Yield (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_seasonal_water_yield.bat" "" "${INVEST_ICON}"

  CreateShortCut "${SMPATH}\RouteDEM (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_routedem.bat" "" "${INVEST_ICON}"
  CreateShortCut "${SMPATH}\DelineateIt (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_delineateit.bat" "" "${INVEST_ICON}"

  CreateShortCut "${SMPATH}\Recreation (${ARCHITECTURE}).lnk" "${INVEST_DATA}\invest_recreation.bat" "" "${INVEST_ICON}"

  ; Write registry keys for convenient uninstallation via add/remove programs.
  ; Inspired by the example at
  ; nsis.sourceforge.net/A_simple_installer_with_start_menu_shortcut_and_uninstaller
  !define REGISTRY_PATH "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_PUBLISHER} ${PRODUCT_NAME} ${PRODUCT_VERSION}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayName"          "${PRODUCT_NAME} ${PRODUCT_VERSION}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "UninstallString"      "${UNINSTALL_PATH}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "QuietUninstallString" "${UNINSTALL_PATH} /S"
  WriteRegStr HKLM "${REGISTRY_PATH}" "InstallLocation"      "$INSTDIR"
  WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayIcon"          "${INVEST_ICON}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "Publisher"            "${PRODUCT_PUBLISHER}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "URLInfoAbout"         "${PRODUCT_WEB_SITE}"
  WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayVersion"       "${PRODUCT_VERSION}"
  WriteRegDWORD HKLM "${REGISTRY_PATH}" "NoModify" 1
  WriteRegDWORD HKLM "${REGISTRY_PATH}" "NoRepair" 1


  ; Actually install the information we want to disk.
  SetOutPath "$INSTDIR"
  File ..\..\LICENSE.txt
  File /nonfatal ..\..\doc\users-guide\build\latex\${PDF_NAME}
  file ..\..\HISTORY.rst

  SetOutPath "$INSTDIR\${INVEST_3_FOLDER}\"
  File /r /x *.hg* /x *.svn* ..\..\${INVEST_3_FOLDER}\*
  ; runmodel.bat is here to help automate testing the UIs.
  File runmodel.bat

;  SetOutPath "$INSTDIR\${INVEST_3_FOLDER_x64}\"
;  File /r /x *.hg* /x *.svn* ..\${INVEST_3_FOLDER_x64}\*

  SetOutPath "$INSTDIR\documentation"
  File /r /x *.hg* /x *.svn* ..\..\doc\users-guide\build\html\*

  ; If the user has provided a custom data zipfile, unzip the data.
  ${If} $LocalDataZipFile != ""
    nsisunz::UnzipToLog $LocalDataZipFile "$INSTDIR"
  ${EndIf}

  ; Write the install log to a text file on disk.
  StrCpy $0 "$INSTDIR\install_log.txt"
  Push $0
  Call DumpLog

SectionEnd

; Only add this section if we're running the installer on Windows 7 or below.
; See InVEST Issue #3515.
; This section is disabled in .onInit if we're running Windows 8 or later.
Section "MSVCRT 2008 Runtime (Recommended)" Sec_VCRedist2008
    File vcredist_x86.exe
    ExecWait "vcredist_x86.exe /q"
    Delete vcredist_x86.exe
SectionEnd

Section "uninstall"
  ; Need to enforce execution level as admin.  See
  ; nsis.sourceforge.net/Shortcuts_removal_fails_on_Windows_Vista
  SetShellVarContext all
  rmdir /r "$SMPROGRAMS\${PACKAGE_NAME}"

  ; Delete the installation directory on disk
  rmdir /r "$INSTDIR"

  ; Delete the entire registry key for this version of RIOS.
  DeleteRegKey HKLM "${REGISTRY_PATH}"
SectionEnd

Var LocalDataZip
Var INSTALLER_DIR

!macro downloadFile RemoteFilepath LocalFilepath
    NSISdl::download "${RemoteFilepath}" ${LocalFilepath}
    Pop $R0 ;Get the status of the file downloaded
    StrCmp $R0 "success" got_it failed
    got_it:
       nsisunz::UnzipToLog ${LocalFilepath} "."
       Delete ${LocalFilepath}
       goto done
    failed:
       MessageBox MB_OK "Download failed: $R0 ${RemoteFilepath}. This might have happened because your Internet connection timed out, or our download server is experiencing problems.  The installation will continue normally, but you'll be missing the ${RemoteFilepath} dataset in your installation.  You can manually download that later by visiting the 'Individual inVEST demo datasets' section of our download page at www.naturalcapitalproject.org."
    done:
       ; Write the install log to a text file on disk.
       StrCpy $0 "$INSTDIR\install_data_${LocalFilepath}_log.txt"
       Push $0
       Call DumpLog
!macroend

!macro downloadData Title Filename AdditionalSize
  Section "${Title}"
    AddSize "${AdditionalSize}"

    ; Check to see if the user defined an 'advanced options' zipfile.
    ; If yes, then we should skip all of this checking, since we only want to use
    ; the data that was in that zip.
    ${If} $LocalDataZipFile != ""
        goto end_of_section
    ${EndIf}

    ; Use a local zipfile if it exists in ./sample_data
    ${GetExePath} $INSTALLER_DIR
    StrCpy $LocalDataZip "$INSTALLER_DIR\sample_data\${Filename}"

;    MessageBox MB_OK "zip: $LocalDataZip"
    IfFileExists "$LocalDataZip" LocalFileExists DownloadFile
    LocalFileExists:
        nsisunz::UnzipToLog "$LocalDataZip" "$INSTDIR"
;        MessageBox MB_OK "found it locally"
       goto done
    DownloadFile:
        ;This is hard coded so that all the download data macros go to the same site
        SetOutPath "$INSTDIR"
        !insertmacro downloadFile "http://data.naturalcapitalproject.org/~dataportal/${DATA_LOCATION}/${Filename}" "${Filename}"
      end_of_section:
      SectionEnd
!macroend

SectionGroup /e "InVEST Datasets" SEC_DATA
  ;here all the numbers indicate the size of the downloads in kilobytes
  ;they were calculated by hand by decompressing all the .zip files and recording
  ;the size by hand.
  SectionGroup "Freshwater Datasets" SEC_FRESHWATER_DATA
    !insertmacro downloadData "Freshwater Base Datasets (optional for freshwater models)" "Freshwater.zip" 4710
    !insertmacro downloadData "Hydropower (optional)" "Hydropower.zip" 100
    !insertmacro downloadData "Seasonal Water Yield: (optional)" "seasonal_water_yield.zip" 500000
  SectionGroupEnd

  SectionGroup "Marine Datasets" SEC_MARINE_DATA
    !insertmacro downloadData "Marine Base Datasets (required for many marine models)" "Marine.zip" 1784696
    !insertmacro downloadData "Aquaculture (optional)" "Aquaculture.zip" 856
    !insertmacro downloadData "Coastal Blue Carbon (optional)" "CoastalBlueCarbon.zip" 856
    !insertmacro downloadData "Coastal Protection (optional)" "CoastalProtection.zip" 117760
    !insertmacro downloadData "Fisheries (optional)" "Fisheries.zip" 784
    !insertmacro downloadData "Habitat Risk Assessment (optional)" "HabitatRiskAssess.zip" 8116
    !insertmacro downloadData "Overlap Analysis (optional)" "OverlapAnalysis.zip" 3692
    !insertmacro downloadData "Scenic Quality (optional)" "ScenicQuality.zip" 9421
    !insertmacro downloadData "Wave Energy (required to run model)" "WaveEnergy.zip" 831620
    !insertmacro downloadData "Wind Energy (required to run model)" "WindEnergy.zip" 4804
    !insertmacro downloadData "Recreation (optional)" "recreation.zip" 24
  SectionGroupEnd

  SectionGroup "Terrestrial Datasets" SEC_TERRESTRIAL_DATA
    !insertmacro downloadData "Crop Production (optional)" "CropProduction.zip" 0
    !insertmacro downloadData "GLOBIO (optional)" "globio.zip" 0
    !insertmacro downloadData "Forest Carbon Edge Effect (required for forest carbon edge model)" "forest_carbon_edge_effect.zip" 8270
    !insertmacro downloadData "Carbon (optional)" "carbon.zip" 728
    !insertmacro downloadData "Terrestrial base datasets (optional for many terrestrial)" "Terrestrial.zip" 587776
    !insertmacro downloadData "Habitat Quality (optional)" "HabitatQuality.zip" 160768
    !insertmacro downloadData "Pollination (optional)" "pollination.zip" 176
    !insertmacro downloadData "Scenario Generator: Rule Based (optional)" "ScenarioGenerator.zip" 0
    !insertmacro downloadData "Scenario Generator: Proximity Based (optional)" "scenario_proximity.zip" 7511
  SectionGroupEnd
SectionGroupEnd

Function .onInit
 ${GetOptions} $CMDLINE "/?" $0
 IfErrors skiphelp showhelp
 showhelp:
     MessageBox MB_OK "InVEST: Integrated Valuation of Ecosystem Services and Tradeoffs$\r$\n\
     $\r$\n\
     For more information about InVEST or the Natural Capital Project, visit our \
     website: http://naturalcapitalproject.org/invest$\r$\n\
     $\r$\n\
     Command-Line Options:$\r$\n\
         /?$\t$\t=$\tDisplay this help and exit$\r$\n\
         /S$\t$\t=$\tSilently install InVEST.$\r$\n\
         /D=$\t$\t=$\tSet the installation directory.$\r$\n\
         /DATAZIP=$\t=$\tUse this sample data zipfile.$\r$\n\
         "
     abort
 skiphelp:

 System::Call 'kernel32::CreateMutexA(i 0, i 0, t "InVEST ${VERSION}") i .r1 ?e'
 Pop $R0

 StrCmp $R0 0 +3
   MessageBox MB_OK|MB_ICONEXCLAMATION "An InVEST ${VERSION} installer is already running."
   Abort

  ${ifNot} ${AtMostWin7}
    ; disable the section if we're not running on Windows 7 or earlier.
    ; This section should not execute for Windows 8 or later.
    SectionGetFlags ${Sec_VCRedist2008} $0
    IntOp $0 $0 & ${SECTION_OFF}
    SectionSetFlags ${Sec_VCRedist2008} $0
    SectionSetText ${Sec_VCRedist2008} ""
  ${endIf}

  ; If the user has defined the /DATAZIP flag, set the 'advanced' option
  ; to the user's defined value.
  ${GetOptions} $CMDLINE "/DATAZIP=" $0
  strcpy $LocalDataZipFile $0
FunctionEnd
