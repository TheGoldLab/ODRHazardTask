function [pdat hfit selectTbl CoeffTable NeuTblAll] =getdatafortable(fileName, monkey, behaviorOnly, UnitUse)

datastruct = getADPODR_dataFromFIRA(fileName, monkey, behaviorOnly);

%% Collect relevant data
% PMF
% Independent variable is "signed cues" -- which is cue location
%   signed by corresponding "LLR for switch"
hazards    = unique(datastruct.ecodes.hazard);
hazards(isnan(hazards)) = [];
numHazards = length(hazards);
sCues      = sign(datastruct.ecodes.llr_for_switch).*abs(datastruct.ecodes.sample_id);
cues       = unique(sCues);
cues(isnan(cues)) = [];
numCues    = length(cues);

% Dependent variable is choice switch re: previous TRUE STATE
%   (note that this is not the same as a choice switch, because 
%   the previous choice might have been in error).
Lswitch    = datastruct.ecodes.choice_switch==1;
Lstay      = datastruct.ecodes.choice_switch==0;

% collect data in pdat
%   dim 1 is LLR
%   dim 2 is pmf, cmf T1/T2 choices
%   dim 3 is hazard
%   dim 4 is for data, ns, sem
pdat  = nans(numCues, 3, numHazards, 3); 
Lc1   = datastruct.ecodes.choice==1;
Lc2   = datastruct.ecodes.choice==2;
Ltask = datastruct.ecodes.task_id>=2;
Lgood = Ltask & datastruct.ecodes.score>=0;
if ~any(Lgood)
    axes(axs(1)); 
    title(sprintf('%s: %d trials', strrep(datastruct.fileName, '_', '-'), sum(Lgood)))
    fprintf('%s: NO TRIALS', datastruct.fileName)
    return
end

    
% For each cue location
for ll = 1:numCues
    Lcue = Lgood & sCues==cues(ll);
    
    % For each hazard
    for hh = 1:numHazards
        Lh = Lcue & datastruct.ecodes.hazard == hazards(hh);
        
        % collect data
        pdat(ll,:,hh,1) = [ ...
            sum(Lh&Lswitch)./sum(Lh&(Lswitch|Lstay)).*100, ...
            nanmean(datastruct.ecodes.RT(Lh&Lc1)), ...
            nanmean(datastruct.ecodes.RT(Lh&Lc2))];
        
        % collect ns
        pdat(ll,:,hh,2) = [ ...
            sum(Lh&(Lswitch|Lstay)), sum(Lh&Lc1), sum(Lh&Lc2)];

        % collect sems
        pdat(ll,2:3,hh,3) = [ ...
            nanse(datastruct.ecodes.RT(Lh&Lc1)), ...
            nanse(datastruct.ecodes.RT(Lh&Lc2))];
    end
end

%% 1. PMF
% First get logistic fit
%   input matrix is h1_bias, h2_bias, llr for switch, choice switch
ldat = cat(2, zeros(sum(Lgood), 2), sCues(Lgood), datastruct.ecodes.choice_switch(Lgood));
for hh = 1:numHazards
    ldat(datastruct.ecodes.hazard(Lgood)==hazards(hh), hh) = 1;
end
ldat = ldat(isfinite(ldat(:,4)),:);
fits = logist_fit(ldat, 'lu1');

hfit= nan(1,2);
for hh = 1:numHazards
    hfit(hh) =  1./(1+exp(-fits(hh)));
end


%%

UseUnitName = str2num(char(UnitUse));
[ a UseUnitInd c] = intersect(datastruct.spikeid,UseUnitName);

sampFreq = 1000;

ub = 1.5;
lb = -1.5;

addpath(genpath('C:\Users\alice\Documents\MATLAB\UsefulAnalysisCode'))

% Get all target configurations
Ts = table2array(datastruct.ecodes(:, {'t1_x', 't1_y', 't2_x', 't2_y'}));
gr = 0.9.*ones(1,3);
tc = [153 216 201]./255;
axlm = 18;

Ltr = datastruct.ecodes.task_id==1 & datastruct.ecodes.score==1;    
% targets    
unTs = unique(Ts(Ltr,1:2), 'rows');


vis = round(nanmean(datastruct.timing.sample_on-datastruct.timing.fp_off)-lb*1000):...
    round(nanmean(datastruct.timing.target_off-datastruct.timing.fp_off)-lb*1000);
mem = round(nanmean(datastruct.timing.target_off-datastruct.timing.fp_off)-lb*1000):...
    round(nanmean(datastruct.timing.fp_off-datastruct.timing.fp_off)-lb*1000);
sac = round(nanmean(datastruct.timing.sac_on-50-datastruct.timing.fp_off)-lb*1000):...
    round(nanmean(datastruct.timing.sac_on+50-datastruct.timing.fp_off)-lb*1000);

rew = round(nanmean(datastruct.timing.sac_on-datastruct.timing.targ_acq)--3000):...
    round(nanmean(datastruct.timing.targ_acq+600-datastruct.timing.targ_acq)--3000);

dec = round(nanmean(datastruct.timing.sample_on-datastruct.timing.fp_off)-lb*1000):...
        round(nanmean(datastruct.timing.fp_off-datastruct.timing.fp_off)-lb*1000);

selectTbl = [];


for u = 1:length(UseUnitName)
    allSpikes = zeros(height(datastruct.timing),(ub-lb).*sampFreq);
    rewallSpikes = zeros(height(datastruct.timing),3701);
    allSpikeTS = [];
    
    SpikeRastmat = nan(height(datastruct.timing),max(max(cellfun(@length,datastruct.spikes))));
    EventRastmat = nan(height(datastruct.timing),length(datastruct.timing{1,:}));
    rewSpikeRastmat = nan(height(datastruct.timing),max(max(cellfun(@length,datastruct.spikes))));
    rewEventRastmat = nan(height(datastruct.timing),length(datastruct.timing{1,:}));
    
    spikeRatevis = nan(1,8);
    spikeRatemem = nan(1,8);
    spikeRatesacc = nan(1,8);
    
    for tr = 1:length(datastruct.ecodes.trial_num)
        fpoff = datastruct.timing.fp_off(tr,:);
        trialTS = datastruct.timing{tr,:};
        trialLockTS = trialTS-fpoff;
        trialSpikeTS = datastruct.spikes{tr,UseUnitInd(u)};
        trialLockSpikeTS = trialSpikeTS-fpoff;
        allSpikeTS = [allSpikeTS; trialLockSpikeTS];
        
        matind = round((trialLockSpikeTS(trialLockSpikeTS>=lb*sampFreq&trialLockSpikeTS<=ub*sampFreq)-lb*sampFreq)+1);%+1 so no zero index
        allSpikes(tr,matind) =1;
        
        SpikeRastmat(tr,1:length(trialLockSpikeTS)) = trialLockSpikeTS;
        EventRastmat(tr,1:length(trialLockTS)) = trialLockTS;
    end
    
    binEdges = 1:100:size(allSpikes,2);
    binnedSpikes = nan(size(allSpikes,1),length(binEdges));
    for b = 1:length(binEdges)
        if b == length(binEdges)
            binnedSpikes(:,b) = sum(allSpikes(:,binEdges(b):end),2);
        else
            binnedSpikes(:,b) = sum(allSpikes(:,binEdges(b):binEdges(b+1)),2);
        end
    end

    allBinnedSpikes{u} = binnedSpikes;
    
    for targ = 1:8
       targInds = find(ismember(Ts(:,1:2)==unTs(targ,:),[1 1],'rows'));
       spikeRatevis(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),vis),2))/(length(vis)/sampFreq);
       spikeRatemem(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),mem),2))/(length(mem)/sampFreq);
       spikeRatesacc(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),sac),2))/(length(sac)/sampFreq);
    %     plot(targ,mean(sum(allSpikes(TrialTypeList{targ},:),2)),'+','Color', targColors(targ,:))
    end
    
    theta = [0:45:315].*pi./180;
    f = @(ParamsVM)errVonMises(theta, spikeRatevis,ParamsVM);
    bound_low = [0 0 0];
    bound_high = [100 359*pi/180 359*pi/180];
    [ParamsVM]=fmincon(f,[0 0 0], [],[],[],[],bound_low, bound_high);
    A_1=ParamsVM(1);
    k_1=ParamsVM(2);
    phi_1=ParamsVM(3);
    theta_ext = [0:359].*pi./180;
    fits_vis_dir = vonMises(theta_ext,A_1,k_1,phi_1);
    
    f = @(ParamsVM)errVonMises(theta, spikeRatemem,ParamsVM);
    bound_low = [0 0 0];
    bound_high = [100 359*pi/180 359*pi/180];
    [ParamsVM]=fmincon(f,[0 0 0], [],[],[],[],bound_low, bound_high);
    A_2=ParamsVM(1);
    k_2=ParamsVM(2);
    phi_2=ParamsVM(3);
    fits_mem_dir = vonMises(theta_ext,A_2,k_2,phi_2);
    
    f = @(ParamsVM)errVonMises(theta, spikeRatesacc,ParamsVM);
    bound_low = [0 0 0];
    bound_high = [100 359*pi/180 359*pi/180];
    [ParamsVM]=fmincon(f,[0 0 0], [],[],[],[],bound_low, bound_high);
    A_3=ParamsVM(1);
    k_3=ParamsVM(2);
    phi_3=ParamsVM(3);
    fits_sac_dir = vonMises(theta_ext,A_3,k_3,phi_3);
    
    dir_vis = round(phi_1*180/pi);
    dir_mem = round(phi_2*180/pi);
    dir_sac = round(phi_3*180/pi);
    unit = UseUnitName(u);
    addselect = table(unit, dir_vis, k_1, dir_mem, k_2, dir_sac, k_3);
    selectTbl = vertcat(selectTbl,addselect);
    
    clear allSpikes allSpikeTS spikeRatevis spikeRatemem spikeRatesacc
end
%%
%Use 'good' trials

TrialN = datastruct.ecodes.trial_num(Lgood);
Switch = Lswitch(Lgood);
Hs = datastruct.ecodes.hazard(Lgood);
CueLoc = datastruct.ecodes.sample_id(Lgood);
CueSw = sCues(Lgood);
prevCorrInd = find(Lgood)-1;
pad = 0;
if sum(prevCorrInd==0)>0
    prevCorrInd(prevCorrInd==0)=[];
    pad = 1;
end
PrevCorr = datastruct.ecodes.score(prevCorrInd);
if pad == 1;
    PrevCorr = [nan; PrevCorr];
end


decbin = find(binEdges<dec(1),1,'last'):find(binEdges>dec(end),1);

CoeffTable =[];
NeuTblAll = [];

for u = 1:length(UseUnitName)
    decSumSpikes = sum(allBinnedSpikes{u}(Lgood,decbin),2);
    decSpikesPerSec = decSumSpikes./length(decbin)./0.1;
    neuTable = table(TrialN,Switch,Hs,CueLoc,CueSw,PrevCorr, decSpikesPerSec);
    
    mdl1 = fitlm(neuTable,'PredictorVars',[1:6],'ResponseVar',"decSpikesPerSec");

    sig = mdl1.Coefficients.pValue(2:end)<0.05;
    coefs = mdl1.Coefficients.Estimate(2:end);
    varnames = mdl1.CoefficientNames(2:end);
    
    neuTable.unit = categorical(UseUnitName(u)*ones(height(neuTable),1));
    NeuTblAll = vertcat(NeuTblAll,neuTable);
    
    tbladd = table(string(varnames'),coefs,sig);
    tbladd.unit = categorical(UseUnitName(u)*ones(height(tbladd),1));
    CoeffTable = vertcat(CoeffTable, tbladd);


    clear decSumSpikes neuTable mdl1
end

%%

end