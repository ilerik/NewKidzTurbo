#!MC 1410
$!VarSet |MFBD| = 'C:\Users\Erik\Dropbox\Science\!Projects\GitHub\NewKidzTurbo'
$!READDATASET  '"STANDARDSYNTAX" "1.0" "FILELIST_CGNSFILES" "1" "dataT0.00000000.cgns" "LoadBCs" "Yes" "AssignStrandIDs" "Yes"'
  DATASETREADER = 'CGNS Loader'
$!ACTIVELINEMAPS -= [1]
$!ACTIVELINEMAPS += [9]
$!VIEW FIT
$!RemoveVar |MFBD|