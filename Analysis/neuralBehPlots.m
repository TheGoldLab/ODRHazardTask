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

dec = round(mean(datastruct.timing.sample_on-datastruct.timing.fp_off)-lb*1000):...
        round(mean(datastruct.timing.fp_off-datastruct.timing.fp_off)-lb*1000);


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
tiledlayout(length(UseUnitName),4);

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
binnedSpikes = nan(size(allSpikes,1),length(binEdges));
for b = 1:length(binEdges)
    if b == length(binEdges)
        binnedSpikes(:,b) = sum(allSpikes(:,binEdges(b):end),2);
    else
        binnedSpikes(:,b) = sum(allSpikes(:,binEdges(b):binEdges(b+1)),2);
    end
end

CueFR = zeros(numCues,length(binEdges));
HFR  = zeros(numHazards,length(binEdges));
CueHFR = zeros(numCues*numHazards,length(binEdges));
normcueHPR = nan(numHazards,numCues);
ch = 1;
decbin = find(binEdges<dec(1),1,'last'):find(binEdges>dec(end),1);
% For each cue location
for ll = 1:numCues
    Lcue = Lgood & sCues==cues(ll);
    CueFR(ll,:) = mean(binnedSpikes(Lcue,:))./0.1;
    
    % For each hazard
    for hh = 1:numHazards
        if ll == 1
            HFR(hh,:) = mean(binnedSpikes(datastruct.ecodes.hazard == hazards(hh),:))./0.1;
        end
        Lh = Lcue & datastruct.ecodes.hazard == hazards(hh);
        CueHFR(ch,:) = mean(binnedSpikes(Lh,:))./0.1;
        ch = ch+1;
        
        
        
        normcueHPR(hh,ll) = sum(sum(binnedSpikes(Lh,decbin)))./...
            sum(sum(binnedSpikes(datastruct.ecodes.hazard == hazards(hh),decbin)));
%         normcueHFR(hh,ll) = sum(sum(allSpikes(Lh,:)))./sum(length(binnedSpikes(Lh,:)));
        
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

C = repmat(linspace(1,0.4,size(CueFR,1)).',1,3);
% axes('ColorOrder',C,'NextPlot','add')
% colormap(gca,'gray')
hold on
for ll = 1:numCues
plot(binEdges+lb*1000,CueFR(ll,:),'Color', C(ll,:))
end
ylabel('FR (sp/s)')
xlabel('Time from fix off (ms)')
title(['FR by Cue - Unit ' num2str(UseUnitName(u))])
axis square
legend(num2str(cues),'Location','BestOutside')

nexttile
plot(binEdges+lb*1000,HFR')
title(['FR by Haz - Unit ' num2str(UseUnitName(u))])
axis square
legend(num2str(hazards),'Location','BestOutside')
nexttile
plot(binEdges+lb*1000,CueHFR')
title(['FR by Haz and Cue - Unit ' num2str(UseUnitName(u))])
axis square
% colormap jet
nexttile
plot(cues,normcueHPR)
xlabel('cue')
ylabel('Percent Response')
title(['Spike Dec Resp to Cue by Block - Unit ' num2str(UseUnitName(u))])
axis square
legend(num2str(hazards),'Location','BestOutside')


end


f1.WindowState = 'maximize';
text(-20,10.6,fileName)

cd('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Figures\SortedSessions')
cd(fileName)
exportgraphics(f1,[fileName '_neuralbeh1.png'],'Resolution',300)
close(f1)

