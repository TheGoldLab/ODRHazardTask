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

figure
plot([0.005 0.5],hfitall','ko-')
xlim([-0.1 0.6])
ylabel('Fit H')
xlabel('Generative H')

figure
subplot(1,3,1)
scatter(selectTblAll.dir_vis,selectTblAll.k_1)
axis square
xlabel('Degrees')
ylabel('k fit')
title('Vis. Select.')
subplot(1,3,2)
scatter(selectTblAll.dir_mem,selectTblAll.k_2)
axis square
xlabel('Degrees')
ylabel('k fit')
title('Mem. Select.')
subplot(1,3,3)
scatter(selectTblAll.dir_sac,selectTblAll.k_3)
axis square
xlabel('Degrees')
ylabel('k fit')
title('Sacc. Select.')

