using MAT
using JSON: JSON

if isdir("FitCheck") == false  mkdir("FitCheck");  end

envVar = JSON.parsefile("param.json");
res_on = get(envVar,"res_on", false);  
const_on = get(envVar,"const_on", false);
vis_on = get(envVar,"vis_on", false);
DirOfBOdotJSON = get(envVar,"DirOfBOdotJSON", "");

DirOfExpData  = get(envVar,"DirOfExpData", "");   
DirOfGridData = get(envVar,"DirOfGridData", "");  
DirOfResData  = get(envVar,"DirOfResData", "");   
DirOfJsonData = get(envVar,"DirOfJsonData", "");

if DirOfExpData[end]  != "/"  DirOfExpData  = DirOfExpData*"/"   end
if DirOfGridData[end] != "/"  DirOfGridData = DirOfGridData*"/"  end
if DirOfResData[end]  != "/"  DirOfResData  = DirOfResData*"/"   end
if DirOfJsonData[end] != "/"  DirOfJsonData = DirOfJsonData*"/"  end

NameOfExpData  = get(envVar,"NameOfExpData", "");
NameOfGridData = get(envVar,"NameOfGridData", "");
NameOfResData  = get(envVar,"NameOfResData", "");

#################################################################################################
#################################################################################################

# =======================================
# Experimental data is loaded here.
# =======================================
AbsExpData  = matopen(DirOfExpData*NameOfExpData);
AbsGridData = matopen(DirOfGridData*NameOfGridData);
AbsResData  = matopen(DirOfResData*NameOfResData);

ExpData02 = read(AbsExpData, "ExpData02");
path02, qGrd02, qItv02, qPtN02, eBot02, eTop02, us02, vs02, ws02, qWid02 = 
  loadGridData(AbsGridData,"02"; type = "spaghettiPlot_v3");

# ! # path02, qGrd02, qItv02, eBot02, eTop02 = loadGridData(AbsGridData,"02"; type = "spaghettiPlot");
# ! # --> this is for old version of spaghetti plot

ExpData03 = read(AbsExpData, "ExpData03");
xBas03, yBas03, zBas03, xGrd03, yGrd03, zSum03, eSum03, HKL03 = loadGridData(AbsGridData,"03");

ExpData04 = read(AbsExpData, "ExpData04");
xBas04, yBas04, zBas04, xGrd04, yGrd04, zSum04, eSum04, HKL04 = loadGridData(AbsGridData,"04");

ExpData05 = read(AbsExpData, "ExpData05");
xBas05, yBas05, zBas05, xGrd05, yGrd05, zSum05, eSum05, HKL05 = loadGridData(AbsGridData,"05");

ExpData06 = read(AbsExpData, "ExpData06");
xBas06, yBas06, zBas06, xGrd06, yGrd06, zSum06, eSum06, HKL06 = loadGridData(AbsGridData,"06");

ExpData07 = read(AbsExpData, "ExpData07");
xBas07, yBas07, zBas07, xGrd07, yGrd07, zSum07, eSum07, HKL07 = loadGridData(AbsGridData,"07");

ExpData08 = read(AbsExpData, "ExpData08");
xBas08, yBas08, zBas08, xGrd08, yGrd08, zSum08, eSum08, HKL08 = loadGridData(AbsGridData,"08");

ExpData09 = read(AbsExpData, "ExpData09");
xBas09, yBas09, zBas09, xGrd09, yGrd09, zSum09, eSum09, HKL09 = loadGridData(AbsGridData,"09");

ExpData10 = read(AbsExpData, "ExpData10");
xBas10, yBas10, zBas10, xGrd10, yGrd10, zSum10, eSum10, HKL10 = loadGridData(AbsGridData,"10");

# =======================================
# Q-grid is generated here.
# Each data has different Q-grid.
# =======================================
qPathData02 = makeSweepPoints(path02,qPtN02,vs02,ws02,qWid02); 
# ? This is for new version of spaghetti plot
# ! # qPathData02 = makeSweepPath(crystl,path02,qGrd02,qItv02);
# ! # --> this is for old version of spaghetti plot

qpts03 = makeGridPoint(xBas03,yBas03,xGrd03,yGrd03,zBas03,zSum03,HKL03);
qpts04 = makeGridPoint(xBas04,yBas04,xGrd04,yGrd04,zBas04,zSum04,HKL04);
qpts05 = makeGridPoint(xBas05,yBas05,xGrd05,yGrd05,zBas05,zSum05,HKL05);
qpts06 = makeGridPoint(xBas06,yBas06,xGrd06,yGrd06,zBas06,zSum06,HKL06);
qpts07 = makeGridPoint(xBas07,yBas07,xGrd07,yGrd07,zBas07,zSum07,HKL07);
qpts08 = makeGridPoint(xBas08,yBas08,xGrd08,yGrd08,zBas08,zSum08,HKL08);
qpts09 = makeGridPoint(xBas09,yBas09,xGrd09,yGrd09,zBas09,zSum09,HKL09);
qpts10 = makeGridPoint(xBas10,yBas10,xGrd10,yGrd10,zBas10,zSum10,HKL10);
