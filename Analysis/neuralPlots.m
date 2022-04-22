
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

    

%For graphing
unitColors = [0.5 0.5 0.5; 0 0 1; 1 0 0; 0 1 0];

f1 = figure;
tiledlayout(length(UseUnitName),3);
f2 = figure;
tiledlayout(length(UseUnitName),3);
f3 = figure;
tiledlayout(length(UseUnitName),2);


for u = 1:length(UseUnitName)
    allSpikes = zeros(height(datastruct.timing),(ub-lb).*sampFreq);
    rewallSpikes = zeros(height(datastruct.timing),3701);
    allSpikeTS = [];
    MGScorrSpikeTS = [];
    MGSerrSpikeTS = [];
    taskcorrSpikeTS = [];
    taskerrSpikeTS = [];
    rewcorrSpikeTS =[];
    rewerrSpikeTS = [];
    SpikeRastmat = nan(height(datastruct.timing),max(max(cellfun(@length,datastruct.spikes))));
    EventRastmat = nan(height(datastruct.timing),length(datastruct.timing{1,:}));
    rewSpikeRastmat = nan(height(datastruct.timing),max(max(cellfun(@length,datastruct.spikes))));
    rewEventRastmat = nan(height(datastruct.timing),length(datastruct.timing{1,:}));
    
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

        targacqrew = datastruct.timing.targ_acq(tr,:);
        rewLockTS = trialTS-targacqrew;
        rewLockSpikeTS = trialSpikeTS-targacqrew;
        if datastruct.ecodes.task_id(tr)>=2
            if datastruct.ecodes.score(tr) == 1
                rewcorrSpikeTS = [rewcorrSpikeTS; rewLockSpikeTS];
            elseif datastruct.ecodes.score(tr) == 0
                rewerrSpikeTS = [rewerrSpikeTS; rewLockSpikeTS];
            end
        end
        
        rewSpikeRastmat(tr,1:length(rewLockSpikeTS)) = rewLockSpikeTS;
        rewEventRastmat(tr,1:length(rewLockTS)) = rewLockTS;
        
        rewmatind = round((rewLockSpikeTS(rewLockSpikeTS>=-3000&rewLockSpikeTS<=700)--3000)+1);%+1 so no zero index
        rewallSpikes(tr,rewmatind) =1;
        
        
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
    
%     rewbin = find(binEdges<rew(1),1,'last'):find(binEdges>rew(end),1);
    
    rewbinEdges = 1:100:size(rewallSpikes,2);
    rewbinnedSpikes = nan(size(rewallSpikes,1),length(rewbinEdges));
    for b = 1:length(rewbinEdges)
        if b == length(rewbinEdges)
            rewbinnedSpikes(:,b) = sum(rewallSpikes(:,rewbinEdges(b):end),2);
        else
            rewbinnedSpikes(:,b) = sum(rewallSpikes(:,rewbinEdges(b):rewbinEdges(b+1)),2);
        end
    end

    rewallBinnedSpikes{u} = rewbinnedSpikes;
    
    rewbin = find(rewbinEdges<rew(1),1,'last'):find(rewbinEdges>rew(end),1);
    
    
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
    
%     nexttile
%     histogram(MGScorrSpikeTS)
%     hold on
%     histogram(MGSerrSpikeTS)
%     xlabel('Time (ms) after fixation point turned off')
%     xlim([-2.5*1000 2.5*1000])
%     title(['MGS PSTH - Unit ' num2str(UseUnitName(u))])
%     axis square
    
    nexttile
    hold on
    plot(linspace(lb*sampFreq,ub*sampFreq,size(binnedSpikes,2)),...
        mean(binnedSpikes(datastruct.ecodes.task_id==1&datastruct.ecodes.score==1,:))./0.1,'k')
    plot(linspace(lb*sampFreq,ub*sampFreq,size(binnedSpikes,2)),...
        mean(binnedSpikes(datastruct.ecodes.task_id==1&datastruct.ecodes.score==-1,:))./0.1,'r')
    xline(0,':')
    xlabel('Time (ms) after fix off')
    ylabel('spikes/s')
    title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
    pbaspect([2,1,1])
    
%     nexttile
%     histogram(taskcorrSpikeTS)
%     hold on
%     histogram(taskerrSpikeTS)
%     xlabel('Time (ms) after fixation point turned off')
%     xlim([-2.5*1000 2.5*1000])
%     title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
%     axis square
    
    nexttile
    hold on
    plot(linspace(lb*sampFreq,ub*sampFreq,size(binnedSpikes,2)),...
        mean(binnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==1,:))./0.1,'k')
    plot(linspace(lb*sampFreq,ub*sampFreq,size(binnedSpikes,2)),...
        mean(binnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==0,:))./0.1,'r')
    xline(0,':')
    xlabel('Time (ms) after fix off')
    ylabel('spikes/s')
    title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
    pbaspect([2,1,1])
    
    
       
    
    figure(f2)

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


    nexttile
    polar([0:45:315].*pi./180,spikeRatevis)
    hold on
    polar(theta_ext,fits_vis_dir)
    text(-1,0,['k = ' num2str(round(k_1,2))])
    text(-2,-A_1*0.4,['dir = ' num2str(round(phi_1*180/pi))])
    title(['Vis. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    polar([0:45:315].*pi./180,spikeRatemem)
    hold on
    polar(theta_ext,fits_mem_dir)
    text(-1,0,['k = ' num2str(round(k_2,2))])
    text(-2,-A_2*0.4,['dir = ' num2str(round(phi_2*180/pi))])
    title(['Mem. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    nexttile
    polar([0:45:315].*pi./180,spikeRatesacc)
    hold on
    polar(theta_ext,fits_sac_dir)
    text(-1,0,['k = ' num2str(round(k_2,2))])
    text(-2,-A_2*0.4,['dir = ' num2str(round(phi_2*180/pi))])
    title(['Sacc. Dir. Select. - Unit ' num2str(UseUnitName(u))])
    axis square
    
    
    figure(f3)
    
    nexttile
    plot(rewEventRastmat,repmat([1:length(datastruct.ecodes.trial_num)]',1,size(rewEventRastmat,2)),'r|', 'MarkerSize', 0.1)
    hold on
    plot(rewSpikeRastmat,repmat([1:length(datastruct.ecodes.trial_num)]',1,size(rewSpikeRastmat,2)),'k|','MarkerSize', 0.1)
    blockChangeInd = find(diff(datastruct.ecodes.task_id)~=0)+1;
    yline(blockChangeInd,'r','LineWidth',1);
    xlim([-3000 700])
    xlabel('Time (ms) after reward expected')
    ylabel('Trial number')
    title(['Raster - Unit ' num2str(UseUnitName(u))])
    axis square
    
%     nexttile
%     histogram(rewcorrSpikeTS)
%     hold on
%     histogram(rewerrSpikeTS)
%     xlabel('Time (ms) after reward expected')
%     xlim([-3*1000 700])
%     title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
%     axis square
    
    nexttile
    hold on
    plot(linspace(-3000,700,size(rewbinnedSpikes,2)),...
        mean(rewbinnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==1,:))./0.1,'k')
    plot(linspace(-3000,700,size(rewbinnedSpikes,2)),...
        mean(rewbinnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==0,:))./0.1,'r')
    xline(0,':')
    xlabel('Time (ms) after reward expected')
    ylabel('spikes/s')
    title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
    pbaspect([3,1,1])
    
%     nexttile
%     hold on
%     plot(linspace(-100*(length(rewbin)-7),600,length(rewbin)),...
%         mean(rewbinnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==1,rewbin))./0.1,'k')
%     plot(linspace(-100*(length(rewbin)-7),600,length(rewbin)),...
%         mean(rewbinnedSpikes(datastruct.ecodes.task_id>1&datastruct.ecodes.score==0,rewbin))./0.1,'r')
%     xline(0,':')
%     xlabel('Time (ms) after reward expected')
%     ylabel('spikes/s')
%     title(['Task PSTH - Unit ' num2str(UseUnitName(u))])
%     axis square
    
   
    
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


f3.WindowState = 'maximize';
sgtitle(fileName)

cd(figLoc)
% mkdir(fileName)
cd(fileName)
exportgraphics(f3,[fileName '_4_rewardresponse.png'],'Resolution',300)
close(f3)
end