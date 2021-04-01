%Get the file
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';

% data_path = '/Users/lab/Desktop/Ci_training/MatFiles';
% figPath = '/Users/lab/Desktop/Ci_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);
f1=figure;
colors = colormap('lines');
f2=figure;
f3=figure;
f4=figure;
f5=figure;

numLookAt = 7;

for f = 1:numLookAt
    currentFile = nfiles-f+1;
    theFile= files(currentFile).name;
    temp=load(theFile);
    data=temp.data;
    
    name = theFile(end-8:end-4);
    name(strfind(name,'_'))='-';
    
    pcorrCue =[];
    
    CueAng = unique(data.ecodes.data(:,38));
    for c = 1:length(CueAng)
       cueInd = data.ecodes.data(:,38)==CueAng(c);
       pcorrCue(c) =  sum(data.ecodes.data(cueInd,41)==1)./(sum(data.ecodes.data(cueInd,41)==1)+sum(data.ecodes.data(cueInd,41)==0));  
    end
    CueAngCon = CueAng;
    if sum(CueAng>270)>0
        CueAngLeftConv = CueAng(CueAng>270)-360;
        CueAngCon(CueAng>270)= CueAngLeftConv;
    end
    CueAngCorr = [CueAngCon pcorrCue'];
    CueAngCorr =sortrows(CueAngCorr);
    figure(f1)
    hold on
    plot(CueAngCorr(:,1),CueAngCorr(:,2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    ylim([0 1])
    pbaspect([3 2 1])
    xlabel('Cue (deg)')
    ylabel('Percent Correct')
    title('P(correct) by cue angle across sessions')
    
    CueType = unique(data.ecodes.data(:,30));
    for c = 1:length(CueType)
       cueInd = data.ecodes.data(:,30)==CueType(c);
       pcorrCueT2(c) =  sum(data.ecodes.data(cueInd,41)==1)./(sum(data.ecodes.data(cueInd,41)==1)+sum(data.ecodes.data(cueInd,41)==0));     
    end
    figure(f2)
    if sum(CueType<205)>0
        hold on
        plot(CueType(CueType<205),pcorrCueT2(CueType<205),'-o','Color',colors(f,:),'DisplayName',...
        name)
    end
    if sum(CueType>=205&CueType<300)>0
        hold on
        plot(CueType(CueType>=205&CueType<300),pcorrCueT2(CueType>=205&CueType<300),'-o','Color',colors(f,:),'DisplayName',...
        name)
    end
    if sum(CueType>=300)>0
        hold on
        plot(CueType(CueType>=300)-80,pcorrCueT2(CueType>=300),'-o','Color',colors(f,:),'DisplayName',...
        name)
    end
    ylim([0 1])
    xlim([190 240])
    pbaspect([7 1 1])
    xlabel('Cue (name)')
    ylabel('Percent Correct')
    title('P(correct) by cue across sessions')
    
    figure(f3)
    corrTrials = data.ecodes.data(~isnan(data.ecodes.data(:,13)),34);
    actTarg = data.ecodes.data(~isnan(data.ecodes.data(:,13)),35);
    actTargPlot = ~isnan(actTarg)&(actTarg==135|actTarg==45);
    subplot(numLookAt,1,f)
    plot(movmean(corrTrials,25),'-','Color',colors(f,:))
    hold on
    plot(actTargPlot,'s')
    pbaspect([5 1 1])
    xlabel('Trials')
    ylabel('Percent Correct')
    title(name)
    
    figure(f4)
    centerInd = find(data.ecodes.data(~isnan(data.ecodes.data(:,13)),38)==0|data.ecodes.data(~isnan(data.ecodes.data(:,13)),38)==180);
    corrTrialsCenter = corrTrials(centerInd);
    actTargPlotCenter=actTargPlot(centerInd);
    subplot(numLookAt,1,f)
    plot(movmean(corrTrialsCenter,25),'-','Color',colors(f,:))
    hold on
    plot(actTargPlotCenter,'s')
    pbaspect([5 1 1])
    xlabel('Trials')
    ylabel('Percent Correct')
    title(['Center Cue Mov. Avg. ',name])
    
    figure(f5)
    finTrialInd = find(~isnan(data.ecodes.data(:,13)));
    sampfreq=data.analog.acquire_rate(end);
    tTargAcq = data.ecodes.data(finTrialInd,13)-data.ecodes.data(finTrialInd,5);
    for tr = 1:length(finTrialInd)
        acqX(tr) = data.analog.data(finTrialInd(tr),2).values(round(tTargAcq(tr)));
        acqY(tr) = data.analog.data(finTrialInd(tr),3).values(round(tTargAcq(tr)));
    end

    choiceDir = acqY>0;

    cueFin =data.ecodes.data(finTrialInd,38);
    actTargFin = data.ecodes.data(finTrialInd,35);
    actTargFinUp = (actTargFin==135|actTargFin==45);

    pChooseUp =[];
    for c = 1:length(CueAng)
       cueInd =[];
        cueInd = cueFin==CueAng(c);
       pChooseUp(c) =  sum(choiceDir(cueInd))./length(choiceDir(cueInd));   
    end
    CueChoose = [CueAngCon pChooseUp'];
    CueChoose = sortrows(CueChoose);
    subplot(1,3,1)
    hold on
    plot(CueChoose(:,1),CueChoose(:,2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    ylim([0 1])
    pbaspect([1 1 1])
    ylabel('pChooseUp')
    xlabel('cue (deg)')
    title('Session')
    pChooseUp =[];
    for c = 1:length(CueAng)
        cueInd =[];
       taskInd = (data.ecodes.data(finTrialInd,29)==3);
        cueInd = intersect(find(cueFin==CueAng(c)),find(taskInd));
       pChooseUp(c) =  sum(choiceDir(cueInd))./length(choiceDir(cueInd));   
    end
    CueChoose = [CueAngCon pChooseUp'];
    CueChoose = sortrows(CueChoose);
    CueChoose(isnan(CueChoose(:,2)),:)=[];
    subplot(1,3,2)
    hold on
    plot(CueChoose(:,1),CueChoose(:,2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    ylim([0 1])
    pbaspect([1 1 1])
    title('Random')
    pChooseUp =[];
    CueChoose=[];
    for c = 1:length(CueAng)
        cueInd =[];
       taskInd = (data.ecodes.data(finTrialInd,29)==2)&(~actTargFinUp);
        cueInd = intersect(find(cueFin==CueAng(c)),find(taskInd));
       pChooseUp(c) =  sum(choiceDir(cueInd))./length(choiceDir(cueInd));   
    end
    CueChoose = [CueAngCon pChooseUp'];
    CueChoose = sortrows(CueChoose);
    CueChoose(isnan(CueChoose(:,2)),:)=[];
    subplot(1,3,3)
    plot(CueChoose(:,1),CueChoose(:,2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    ylim([0 1])
    pbaspect([1 1 1])
    hold on
    pChooseUp =[];
    CueChoose =[];
    for c = 1:length(CueAng)
        cueInd =[];
       taskInd = (data.ecodes.data(finTrialInd,29)==2)&(actTargFinUp);
        cueInd = intersect(find(cueFin==CueAng(c)),find(taskInd));
       pChooseUp(c) =  sum(choiceDir(cueInd))./length(choiceDir(cueInd));   
    end
    CueChoose = [CueAngCon pChooseUp'];
    CueChoose = sortrows(CueChoose);
    CueChoose(isnan(CueChoose(:,2)),:)=[];
    plot(CueChoose(:,1),CueChoose(:,2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    title('No change')
end

figure(f1)
legend TOGGLE
legend('location', 'Southeast')

figure(f2)
legend TOGGLE
pbaspect([5 2 1])
legend('location', 'southwest')

figure(f5)
legend TOGGLE
legend('location', 'BestOutside')