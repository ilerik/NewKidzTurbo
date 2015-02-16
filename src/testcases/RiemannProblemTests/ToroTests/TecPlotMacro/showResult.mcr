#!MC 1410
$!VarSet |MFBD| = 'C:\Users\Erik\Dropbox\Science\!Projects\GitHub\NewKidzTurbo'
$!READDATASET  '"STANDARDSYNTAX" "1.0" "FILELIST_CGNSFILES" "1" "result.cgns" "LoadBCs" "Yes" "AssignStrandIDs" "Yes" "LoaderVersion" "V3" "CgnsLibraryVersion" "3.1.4"'
  DATASETREADER = 'CGNS Loader'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  ASSIGNSTRANDIDS = NO
  INITIALPLOTTYPE = XYLINE
  INITIALPLOTFIRSTZONEONLY = NO
  ADDZONESTOEXISTINGSTRANDS = NO
$!ACTIVELINEMAPS -= [1-17]
$!DELETELINEMAPS  [1-17]
$!CREATELINEMAP 
$!LINEMAP [1]  NAME = '&ZN&'
$!ACTIVELINEMAPS += [1]
$!CREATELINEMAP 
$!LINEMAP [2]  NAME = '&ZN&'
$!LINEMAP [2]  ASSIGN{ZONE = 2}
$!ACTIVELINEMAPS += [2]
$!CREATELINEMAP 
$!LINEMAP [3]  NAME = '&ZN&'
$!LINEMAP [3]  ASSIGN{ZONE = 3}
$!ACTIVELINEMAPS += [3]
$!CREATELINEMAP 
$!LINEMAP [4]  NAME = '&ZN&'
$!LINEMAP [4]  ASSIGN{ZONE = 4}
$!ACTIVELINEMAPS += [4]
$!VIEW FIT
$!CREATELINEMAP 
$!LINEMAP [5]  NAME = '&ZN&'
$!LINEMAP [5]  ASSIGN{YAXISVAR = 10}
$!CREATELINEMAP 
$!LINEMAP [6]  NAME = '&ZN&'
$!LINEMAP [6]  ASSIGN{YAXISVAR = 10}
$!LINEMAP [6]  ASSIGN{ZONE = 2}
$!CREATELINEMAP 
$!LINEMAP [7]  NAME = '&ZN&'
$!LINEMAP [7]  ASSIGN{YAXISVAR = 10}
$!LINEMAP [7]  ASSIGN{ZONE = 3}
$!CREATELINEMAP 
$!LINEMAP [8]  NAME = '&ZN&'
$!LINEMAP [8]  ASSIGN{YAXISVAR = 10}
$!LINEMAP [8]  ASSIGN{ZONE = 4}
$!PICK ADDATPOSITION
  X = 4.352248394
  Y = 2.48554603854
  CONSIDERSTYLE = YES
$!VIEW FIT
$!CREATELINEMAP 
$!LINEMAP [9]  NAME = '&ZN&'
$!LINEMAP [9]  ASSIGN{YAXISVAR = 7}
$!CREATELINEMAP 
$!LINEMAP [10]  NAME = '&ZN&'
$!LINEMAP [10]  ASSIGN{YAXISVAR = 7}
$!LINEMAP [10]  ASSIGN{ZONE = 2}
$!CREATELINEMAP 
$!LINEMAP [11]  NAME = '&ZN&'
$!LINEMAP [11]  ASSIGN{YAXISVAR = 7}
$!LINEMAP [11]  ASSIGN{ZONE = 3}
$!CREATELINEMAP 
$!LINEMAP [12]  NAME = '&ZN&'
$!LINEMAP [12]  ASSIGN{YAXISVAR = 7}
$!LINEMAP [12]  ASSIGN{ZONE = 4}
$!PICK ADDATPOSITION
  X = 4.56638115632
  Y = 2.20289079229
  CONSIDERSTYLE = YES
$!VIEW FIT

$!LINEPLOTLAYERS SHOWSYMBOLS = YES
$!PICK ADDATPOSITION
  X = 2.45074946467
  Y = 3.23072805139
  CONSIDERSTYLE = YES
$!VIEW FIT
$!ACTIVELINEMAPS -= [4]
$!PICK ADDATPOSITION
  X = 3.04175588865
  Y = 2.21145610278
  CONSIDERSTYLE = YES
$!VIEW FIT
$!ACTIVELINEMAPS += [4]
$!LINEMAP [4]  SYMBOLS{SHOW = NO}
$!PICK ADDATPOSITION
  X = 0.369379014989
  Y = 2.00588865096
  CONSIDERSTYLE = YES
$!VIEW FIT
$!GLOBALLINEPLOT LEGEND{SHOW = YES}
$!PICK ADDATPOSITION
  X = 1.32869379015
  Y = 3.05942184154
  CONSIDERSTYLE = YES
$!VIEW FIT
$!LINEMAP [1]  LINES{SHOW = NO}
$!LINEMAP [2]  LINES{SHOW = NO}
$!LINEMAP [3]  LINES{SHOW = NO}
$!PICK ADDATPOSITION
  X = 0.8147751606
  Y = 2.94807280514
  CONSIDERSTYLE = YES
$!VIEW FIT
$!LINEMAP [1]  SYMBOLS{SIZE = 1}
$!LINEMAP [2]  SYMBOLS{SIZE = 1}
$!LINEMAP [3]  SYMBOLS{SIZE = 1}
$!PICK ADDATPOSITION
  X = 1.86830835118
  Y = 3.41059957173
  CONSIDERSTYLE = YES
$!VIEW FIT
$!LINEMAP [2]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = DEL}}
$!LINEMAP [3]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = CIRCLE}}
$!PICK ADDATPOSITION
  X = 1.73982869379
  Y = 3.4534261242
  CONSIDERSTYLE = YES
$!VIEW FIT

$!RemoveVar |MFBD|