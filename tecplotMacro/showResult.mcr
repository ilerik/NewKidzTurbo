#!MC 1410
$!VarSet |MFBD| = 'C:\Program Files\Tecplot\Tecplot 360 EX 2015 R1'
$!READDATASET  '"STANDARDSYNTAX" "1.0" "FILELIST_CGNSFILES" "1" "C:\Users\Erik\Dropbox\Science\!Projects\GitHub\NewKidzTurbo\result.cgns" "LoadBCs" "Yes" "AssignStrandIDs" "Yes" "LoaderVersion" "V3" "CgnsLibraryVersion" "3.1.4"'
  DATASETREADER = 'CGNS Loader'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  ASSIGNSTRANDIDS = NO
  INITIALPLOTTYPE = XYLINE
  INITIALPLOTFIRSTZONEONLY = NO
  ADDZONESTOEXISTINGSTRANDS = NO
$!PICK ADDATPOSITION
  X = -0.487152034261
  Y = 1.62901498929
  CONSIDERSTYLE = YES
$!VIEW FIT
$!RemoveVar |MFBD|
