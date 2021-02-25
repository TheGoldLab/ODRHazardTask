%Get the file
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);
f1=figure;
colors = colormap('lines');
f2=figure;


for f = 1:7
    currentFile = nfiles-f+1;
    theFile= files(currentFile).name;
    temp=load(theFile);
    data=temp.data;
    
    name = theFile(end-8:end-4);
    name(strfind(name,'_'))='-';
    
    CueType = unique(data.ecodes.data(:,30));
    for c = 1:length(CueType)
       cueInd = data.ecodes.data(:,30)==CueType(c);
       pcorrCueT2(c) =  sum(data.ecodes.data(cueInd,34))./length(data.ecodes.data(cueInd,34));   
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
    ylim([0 1])
    xlim([200 209])
    xlabel('Cue (name)')
    ylabel('Percent Correct')
    title('P(correct) by cue across sessions')
    legend TOGGLE

    
    figure(f2)
    corrTrials = data.ecodes.data(:,34);
    actTarg = data.ecodes.data(:,35);
    actTargPlot = actTarg==135&~isnan(actTarg);
    subplot(7,1,f)
    plot(movmean(corrTrials,25),'-','Color',colors(f,:))
    hold on
    plot(actTargPlot,'s')
    pbaspect([7 1 1])
    xlabel('Trials')
    ylabel('Percent Correct')
    title(name)
end