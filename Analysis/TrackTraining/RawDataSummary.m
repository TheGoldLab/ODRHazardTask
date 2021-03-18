%Get the file
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);
f1=figure;
colors = colormap('lines');
f2=figure;
f3=figure;
f4=figure;


for f = 1:3
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
       pcorrCue(c) =  sum(data.ecodes.data(cueInd,34))./length(data.ecodes.data(cueInd,34));   
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
    xlabel('Cue (deg)')
    ylabel('Percent Correct')
    title('P(correct) by cue angle across sessions')
    
    CueType = unique(data.ecodes.data(:,30));
    for c = 1:length(CueType)
       cueInd = data.ecodes.data(:,30)==CueType(c);
       pcorrCueT2(c) =  sum(data.ecodes.data(cueInd,34))./sum(~isnan(data.ecodes.data(cueInd,13)));   
    end
    figure(f2)
    hold on
    plot(CueType(1:2),pcorrCueT2(1:2),'-o','Color',colors(f,:),'DisplayName',...
        name)
    if length(CueType)>2
        hold on
        plot(CueType(3:4),pcorrCueT2(3:4),'-o','Color',colors(f,:),'DisplayName',...
        name)
    end
    if length(CueType)>4
        hold on
        plot(CueType(5:6),pcorrCueT2(3:4),'-o','Color',colors(f,:),'DisplayName',...
        name)
    end
    ylim([0 1])
    xlim([200 309])
    xlabel('Cue (name)')
    ylabel('Percent Correct')
    title('P(correct) by cue across sessions')
    
    figure(f3)
    corrTrials = data.ecodes.data(~isnan(data.ecodes.data(:,13)),34);
    actTarg = data.ecodes.data(~isnan(data.ecodes.data(:,13)),35);
    actTargPlot = ~isnan(actTarg)&(actTarg==135|actTarg==45);
    subplot(7,1,f)
    plot(movmean(corrTrials,25),'-','Color',colors(f,:))
    hold on
    plot(actTargPlot,'s')
    pbaspect([7 1 1])
    xlabel('Trials')
    ylabel('Percent Correct')
    title(name)
    
    figure(f4)
    centerInd = find(data.ecodes.data(~isnan(data.ecodes.data(:,13)),38)==0|data.ecodes.data(~isnan(data.ecodes.data(:,13)),38)==180);
    corrTrialsCenter = corrTrials(centerInd);
    actTargPlotCenter=actTargPlot(centerInd);
    subplot(7,1,f)
    plot(movmean(corrTrialsCenter,25),'-','Color',colors(f,:))
    hold on
    plot(actTargPlotCenter,'s')
    pbaspect([7 1 1])
    xlabel('Trials')
    ylabel('Percent Correct')
    title(['Center Cue Mov. Avg. ',name])
end

figure(f1)
legend TOGGLE
legend('location', 'Best')

figure(f2)
legend TOGGLE
legend('location', 'Best')