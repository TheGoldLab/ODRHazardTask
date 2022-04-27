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

for f = [2:18,20:length(filesUse)]%1:length(monkeyInd)
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
% axis square
xlabel('Degrees')
ylabel('k fit')
title('Vis. Select.')
subplot(1,3,2)
scatter(selectTblAll.dir_mem,selectTblAll.k_2)
% axis square
xlabel('Degrees')
ylabel('k fit')
title('Mem. Select.')
subplot(1,3,3)
scatter(selectTblAll.dir_sac,selectTblAll.k_3)
% axis square
xlabel('Degrees')
ylabel('k fit')
title('Sacc. Select.')

PredVars = unique(CoeffTableAll.Var1);
figure
hold on
for c = 1:length(PredVars)
    rowind = contains(CoeffTableAll.Var1,PredVars(c));
    plot(c,CoeffTableAll.coefs(rowind&CoeffTableAll.sig),'k*')
    plot(c,CoeffTableAll.coefs(rowind&~CoeffTableAll.sig),'o','Color',[0.5 0.5 0.5])
end
set(gca,'XTickLabel',PredVars)
ylabel('Estimate')
title('Each Unit LM Coeffs')

mdl1 = fitlm(NeuTblAllSess,'PredictorVars',[1:6],'ResponseVar',"decSpikesPerSec")

figure
subplot(2,1,1)
plot(mdl1)
subplot(2,1,2)
sig = mdl1.Coefficients.pValue(2:end)<0.05;
coefs = mdl1.Coefficients.Estimate(2:end);
[rs,c]=find(sig);
[rns,c]=find(~sig);

hold on
plot(rs,coefs(find(sig)), '*')
plot(rns,coefs(find(~sig)), 'o')
yline(0,'--')
xticks([1:length(sig)])
set(gca,'XTickLabel',mdl1.CoefficientNames(2:end)) 
xtickangle(45)
ylabel('Coefficients Value')
title('Linear Model Coefficients')
sgtitle('Grouped Linear Model')

