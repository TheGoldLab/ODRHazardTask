
addpath(genpath('C:\Users\alice\Documents\Projects\ODRHazardTask\Analysis\JoshProcessing'));

datastruct = getADPODR_dataFromFIRA(fileName, 'MM', 0)


%datastruct.ecodes.Properties.VariableNames

UseUnitInd = find(mod(datastruct.spikeid,10)~=0);
UseUnitName = datastruct.spikeid(UseUnitInd);

sampFreq = 1000;

ub = 1;
lb = -1.5;
sigma = 0.035.*sampFreq;
dtau = linspace(-250,250,0.25.*sampFreq);
w = sampFreq.*normpdf(dtau,0,sigma);

% Get all target configurations
    Ts = table2array(datastruct.ecodes(:, {'t1_x', 't1_y', 't2_x', 't2_y'}));
    gr = 0.9.*ones(1,3);
    tc = [153 216 201]./255;
    axlm = 18;
    
    Ltr = datastruct.ecodes.task_id==1 & datastruct.ecodes.score==1;    
    % targets    
    unTs = unique(Ts(Ltr,1:2), 'rows');
    
    
    vis = round(mean(datastruct.timing.sample_on-datastruct.timing.fp_off)-lb*1000):...
        round(mean(datastruct.timing.target_off-datastruct.timing.fp_off)-lb*1000);
    mem = mean(datastruct.timing.fp_off-datastruct.timing.target_off);
    mot = mean(datastruct.timing.all_off-datastruct.timing.sac_on);

    

%For graphing
unitColors = [0.5 0.5 0.5; 0 0 1; 1 0 0; 0 1 0];

f1 = figure
tiledlayout(length(UseUnitName),5);

for u = 1:length(UseUnitName)
    allSpikes = zeros(height(datastruct.timing),(ub-lb).*sampFreq);
    allSpikeTS = [];
    %spikeRate = zeros(1,length(targ));
    
    
    nexttile
    for tr = 1:length(datastruct.ecodes.trial_num)
        fpoff = datastruct.timing.fp_off(tr,:);
        trialTS = datastruct.timing{tr,:};
        trialLockTS = trialTS-fpoff;
        trialSpikeTS = datastruct.spikes{tr,UseUnitInd(u)};
        trialLockSpikeTS = trialSpikeTS-fpoff;
        allSpikeTS = [allSpikeTS; trialLockSpikeTS];
        matind = round((trialLockSpikeTS(trialLockSpikeTS>=lb*sampFreq&trialLockSpikeTS<=ub*sampFreq)-lb*sampFreq)+1);
        allSpikes(tr,matind) =1;

        eventBar = [0 1]+tr-1;
        plot(repmat(trialLockTS,2,1),repmat(eventBar,length(trialLockTS),1)','r')
        hold on
        spikeBar = [0 1]+tr-1;
        plot(repmat(trialLockSpikeTS',2,1),repmat(spikeBar,numel(trialLockSpikeTS),1)','k')
        xlabel('Time after fixation point turned off')
    end
    xlim([-2.5*1000 2.5*1000])
    title(['Raster - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    histogram(allSpikeTS)
    xlabel('Time after fixation point turned off')
    xlim([-2.5*1000 2.5*1000])
    title(['PSTH - Unit ' num2str(UseUnitName(u))])
    axis square
    

    for targ = 1:8
       targInds = find(ismember(Ts(:,1:2)==unTs(targ,:),[1 1],'rows'))
       spikeRatevis(targ) = mean(sum(allSpikes(Ltr&targInds,(length(allSpikes)/5-0.2*sampFreq):(length(allSpikes)/5-0.095*sampFreq)),2));
       spikeRatemem(targ) = mean(sum(allSpikes(Ltr&targInds,(length(allSpikes)/5):(length(allSpikes)/5*3)),2));
       spikeRatesacc(targ) = mean(sum(allSpikes(Ltr&targInds,(length(allSpikes)/5*3):(length(allSpikes)/5*3+0.3*sampFreq)),2));
    %     plot(targ,mean(sum(allSpikes(TrialTypeList{targ},:),2)),'+','Color', targColors(targ,:))
    end


    nexttile
    polar([0:15:345].*pi./180,spikeRatevis)
    title(['Vis. Dir. Select. - Unit ' num2str(unitNums(u))])
    axis square
    
    nexttile
    polar([0:15:345].*pi./180,spikeRatemem)
    title(['Mem. Dir. Select. - Unit ' num2str(unitNums(u))])
    axis square
    
    nexttile
    polar([0:15:345].*pi./180,spikeRatesacc)
    title(['Sacc. Dir. Select. - Unit ' num2str(unitNums(u))])
    axis square

    
end