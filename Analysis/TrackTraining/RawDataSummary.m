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


for f = 1:4
    currentFile = nfiles-f+1;
    theFile= files(currentFile).name;
    temp=load(theFile);
    data=temp.data;
    
    name = theFile(end-8:end-4);
    name(strfind(name,'_'))='-';
    
    CueType = unique(data.ecodes.data(:,30));
    for c = 1:length(CueType)
       cueInd = data.ecodes.data(:,30)==CueType(c);
       pcorrCueT2(c) =  sum(data.ecodes.data(cueInd,34))./sum(~isnan(data.ecodes.data(cueInd,13)));   
    end
    figure(f1)
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
    
    figure(f2)
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
    
    figure(f3)
    centerInd = find(data.ecodes.data(~isnan(data.ecodes.data(:,13)),38)==180);
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