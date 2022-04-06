addpath(genpath('C:\Users\alice\Documents\Projects\ODRHazardTask\Analysis\JoshProcessing'));

datastruct = getADPODR_dataFromFIRA(fileName, 'MM', 0)


%datastruct.ecodes.Properties.VariableNames

UseUnitInd = find(mod(datastruct.spikeid,10)~=0);
UseUnitName = datastruct.spikeid(UseUnitInd);

sampFreq = 1000;

ub = 1.5;
lb = -1.5;
sigma = 0.035.*sampFreq;
dtau = linspace(-250,250,0.25.*sampFreq);
w = sampFreq.*normpdf(dtau,0,sigma);


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
    


%%
f1 = figure
tiledlayout(length(UseUnitName),2);

for u = 1:length(UseUnitName)
    allSpikes = zeros(height(datastruct.timing),(ub-lb).*sampFreq);
    allSpikeTS = [];

nexttile
    for tr = 1:length(datastruct.ecodes.trial_num)
        fpoff = datastruct.timing.fp_off(tr,:);
        trialTS = datastruct.timing{tr,:};
        trialLockTS = trialTS-fpoff;
        trialSpikeTS = datastruct.spikes{tr,UseUnitInd(u)};
        trialLockSpikeTS = trialSpikeTS-fpoff;
        matind = round((trialLockSpikeTS(trialLockSpikeTS>=lb*sampFreq&trialLockSpikeTS<=ub*sampFreq)-lb*sampFreq)+1);%+1 so no zero index
        allSpikes(tr,matind) =1;
        
    end

    binEdges = 1:100:size(allSpikes,2);
%     [r, c] = find(allSpikes)
%     binnedSpikes = histogram(c,bins)
binnedSpikes = nan(size(allSpikes,1),length(binEdges)-1);
for b = 1:length(binEdges)-1
    binnedSpikes(:,b) = sum(allSpikes(:,binEdges(b):binEdges(b+1)),2);
end

% For each cue location
for ll = 1:numCues
    Lcue = Lgood & sCues==cues(ll);
    
    % For each hazard
    for hh = 1:numHazards
        Lh = Lcue & datastruct.ecodes.hazard == hazards(hh);
        
%         % collect data
%         pdat(ll,:,hh,1) = [ ...
%             sum(Lh&Lswitch)./sum(Lh&(Lswitch|Lstay)).*100, ...
%             nanmean(data.ecodes.RT(Lh&Lc1)), ...
%             nanmean(data.ecodes.RT(Lh&Lc2))];
%         
%         % collect ns
%         pdat(ll,:,hh,2) = [ ...
%             sum(Lh&(Lswitch|Lstay)), sum(Lh&Lc1), sum(Lh&Lc2)];
% 
%         % collect sems
%         pdat(ll,2:3,hh,3) = [ ...
%             nanse(data.ecodes.RT(Lh&Lc1)), ...
%             nanse(data.ecodes.RT(Lh&Lc2))];
    end
end

end