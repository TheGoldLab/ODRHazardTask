%Wrapper for AODR analysis pipeline

%Add appropriate directories
addAODRpathsADpc
AnalysisLoc = pwd;
addpath(genpath(AnalysisLoc));

cd('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Sorted')
MasterList = readtable('SortedUnits_MasterList');

Monkey = 'MM';
monkeyInd = contains(MasterList.Filename,Monkey);

filesUse = MasterList.Filename(monkeyInd);
UnitsUse = MasterList.Unit(monkeyInd);

hfitall = [];
selectTblAll = [];
CoeffTableAll = [];
NeuTblAllSess = [];

for f = [3:18,20:length(filesUse)]%1:length(monkeyInd)
    fileName = filesUse{f}(1:end-4);
    UnitUse = strsplit(UnitsUse{f}, ', ');
    [pdat, hfit, selectTbl, CoeffTable, NeuTblAll] = getdatafortable(fileName, Monkey, 0, UnitUse);
    
    selectTbl.session = repmat(fileName(1:13),length(UnitUse),1);
    CoeffTable.session = repmat(fileName(1:13),height(CoeffTable),1);
    NeuTblAll.session = repmat(fileName(1:13),height(NeuTblAll),1);
    
    hfitall = vertcat(hfitall, hfit);
    selectTblAll = vertcat(selectTblAll, selectTbl);
    CoeffTableAll = vertcat(CoeffTableAll, CoeffTable);
    NeuTblAllSess = vertcat(NeuTblAllSess, NeuTblAll);
    
end