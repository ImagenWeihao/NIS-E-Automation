[Setup]
AppName=SpheroidPA
AppVersion=1.0
AppPublisher=Imagen BioWorks
DefaultDirName={autopf}\SpheroidPA
DefaultGroupName=SpheroidPA
OutputDir=installer_out
OutputBaseFilename=SpheroidPA_Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern

[Files]
Source: "dist\SpheroidPA\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\SpheroidPA"; Filename: "{app}\SpheroidPA.exe"
Name: "{commondesktop}\SpheroidPA"; Filename: "{app}\SpheroidPA.exe"

[Run]
Filename: "{app}\SpheroidPA.exe"; Description: "Launch SpheroidPA"; Flags: nowait postinstall skipifsilent
