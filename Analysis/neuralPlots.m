
function neuralPlots(fileName, monkey, behaviorOnly, figLoc)

datastruct = getADPODR_dataFromFIRA(fileName, monkey, behaviorOnly);


%datastruct.ecodes.Properties.VariableNames

UseUnitInd = find(mod(datastruct.spikeid,10)~=0);
UseUnitName = datastruct.spikeid(UseUnitInd);

sampFreq = 1000;

ub = 1.5;
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
    
    
    vis = round(nanmean(datastruct.timing.sample_on-datastruct.timing.fp_off)-lb*1000):...
        round(nanmean(datastruct.timing.target_off-datastruct.timing.fp_off)-lb*1000);
    mem = round(nanmean(datastruct.timing.target_off-datastruct.timing.fp_off)-lb*1000):...
        round(nanmean(datastruct.timing.fp_off-datastruct.timing.fp_off)-lb*1000);
    sac = round(nanmean(datastruct.timing.sac_on-50-datastruct.timing.fp_off)-lb*1000):...
        round(nanmean(datastruct.timing.sac_on+50-datastruct.timing.fp_off)-lb*1000);

    

%For graphing
unitColors = [0.5 0.5 0.5; 0 0 1; 1 0 0; 0 1 0];

f1 = figure;
tiledlayout(length(UseUnitName),3);
f2 = figure;
tiledlayout(length(UseUnitName),3);

for u = 1:length(UseUnitName)
    allSpikes = zeros(height(datastruct.timing),(ub-lb).*sampFreq);
    allSpikeTS = [];
    MGScorrSpikeTS = [];
    MGSerrSpikeTS = [];
    taskcorrSpikeTS = [];
    taskerrSpikeTS = [];
    SpikeRastmat = nan(height(datastruct.timing),max(max(cellfun(@length,datastruct.spikes))));
    EventRastmat = nan(height(datastruct.timing),length(datastruct.timing{1,:}));
    
    spikeRatevis = nan(1,8);
    spikeRatemem = nan(1,8);
    spikeRatesacc = nan(1,8);
    figure(f1)
    nexttile
    for tr = 1:length(datastruct.ecodes.trial_num)
        fpoff = datastruct.timing.fp_off(tr,:);
        trialTS = datastruct.timing{tr,:};
        trialLockTS = trialTS-fpoff;
        trialSpikeTS = datastruct.spikes{tr,UseUnitInd(u)};
        trialLockSpikeTS = trialSpikeTS-fpoff;
        allSpikeTS = [allSpikeTS; trialLockSpikeTS];
        if datastruct.ecodes.task_id(tr)==1
            if datastruct.ecodes.score(tr) == 1
                MGScorrSpikeTS = [MGScorrSpikeTS; trialLockSpikeTS];
            elseif datastruct.ecodes.score(tr) == -1
                MGSerrSpikeTS = [MGSerrSpikeTS; trialLockSpikeTS];
            end
        elseif datastruct.ecodes.task_id(tr)>=2
            if datastruct.ecodes.score(tr) == 1
                taskcorrSpikeTS = [taskcorrSpikeTS; trialLockSpikeTS];
            elseif datastruct.ecodes.score(tr) == 0
                taskerrSpikeTS = [taskerrSpikeTS; trialLockSpikeTS];
            end
        end
        
        matind = round((trialLockSpikeTS(trialLockSpikeTS>=lb*sampFreq&trialLockSpikeTS<=ub*sampFreq)-lb*sampFreq)+1);%+1 so no zero index
        allSpikes(tr,matind) =1;
        
        SpikeRastmat(tr,1:length(trialLockSpikeTS)) = trialLockSpikeTS;
        EventRastmat(tr,1:length(trialLockTS)) = trialLockTS;

%         eventBar = [0 1]+tr-1;
%         plot(repmat(trialLockTS,2,1),repmat(eventBar,length(trialLockTS),1)','r')
%         hold on
%         spikeBar = [0 1]+tr-1;
%         plot(repmat(trialLockSpikeTS',2,1),repmat(spikeBar,numel(trialLockSpikeTS),1)','k')
%         xlabel('Time after fixation point turned off')
    end
    plot(EventRastmat,repmat([1:length(datastruct.ecodes.trial_num)]',1,size(EventRastmat,2)),'r|', 'MarkerSize', 0.1)
    hold on
    plot(SpikeRastmat,repmat([1:length(datastruct.ecodes.trial_num)]',1,size(SpikeRastmat,2)),'k|','MarkerSize', 0.1)
    blockChangeInd = find(diff(datastruct.ecodes.task_id)~=0)+1;
    yline(blockChangeInd,'r','LineWidth',1);
    xlim([-2.5*1000 2.5*1000])
    xlabel('Time (ms) after fixation point off')
    ylabel('Trial number')
    title(['Raster - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    histogram(MGScorrSpikeTS)
    hold on
    histogram(MGSerrSpikeTS)
    xlabel('Time (ms) after fixation point turned off')
    xlim([-2.5*1000 2.5*1000])
    title(['MGS PSTH - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    histogram(taskcorrSpikeTS)
    hold on
    histogram(taskerrSpikeTS)
    xlabel('Time (ms) after fixation point turned off')
    xlim([-2.5*1000 2.5*1000])
    title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
    axis square
    
    
    figure(f2)

    for targ = 1:8
       targInds = find(ismember(Ts(:,1:2)==unTs(targ,:),[1 1],'rows'));
       spikeRatevis(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),vis),2));
       spikeRatemem(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),mem),2));
       spikeRatesacc(targ) = mean(sum(allSpikes(intersect(find(Ltr),targInds),sac),2));
    %     plot(targ,mean(sum(allSpikes(TrialTypeList{targ},:),2)),'+','Color', targColors(targ,:))
    end


    nexttile
    polar([0:45:315].*pi./180,spikeRatevis)
    title(['Vis. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    polar([0:45:315].*pi./180,spikeRatemem)
    title(['Mem. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    polar([0:45:315].*pi./180,spikeRatesacc)
    title(['Sacc. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    
    clear allSpikes allSpikeTS spikeRatevis spikeRatemem spikeRatesacc
end

f1.WindowState = 'maximize';
sgtitle(fileName)

cd(figLoc)
% mkdir(fileName)
cd(fileName)
exportgraphics(f1,[fileName '_2_neuraloverview.png'],'Resolution',300)
close(f1)

f2.WindowState = 'maximize';
sgtitle(fileName)

cd(figLoc)
% mkdir(fileName)
cd(fileName)
exportgraphics(f2,[fileName '_3_neuralselectivity.png'],'Resolution',300)
close(f2)

end