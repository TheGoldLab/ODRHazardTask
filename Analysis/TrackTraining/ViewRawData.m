%Get the file
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = 7;

theFile= files(currentFile).name;
[~,b,~]=fileparts(theFile);
savePath = b;

temp=load(theFile);
data=temp.data;

%make a save path to  new folder
cd(figPath)

newFileStr = input('Is this a new file?(y/n):','s');
if strcmpi(newFileStr,'y')
    newMatFile
end

b(strfind(b,'_'))='-';

cd(savePath)

% Make some plots

CueAng = unique(data.ecodes.data(:,38));
for c = 1:length(CueAng)
   cueInd = data.ecodes.data(:,38)==CueAng(c);
   pcorrCue(c) =  sum(data.ecodes.data(cueInd,34))./length(data.ecodes.data(cueInd,34));   
end
figure
plot(CueAng,pcorrCue,'-o')
ylim([0 1])
xlabel('Cue (deg)')
ylabel('Percent Correct')
title(b)
% saveas(gcf,['CorrectbyCue',b],'png')
exportgraphics(gcf,['CorrectbyCue',b,'.png'],'Resolution',300)

CueType = unique(data.ecodes.data(:,30));
for c = 1:length(CueType)
   cueInd = data.ecodes.data(:,30)==CueType(c);
   pcorrCueT2(c) =  sum(data.ecodes.data(cueInd,34))./length(data.ecodes.data(cueInd,34));   
end
figure
plot(CueType(1:2),pcorrCueT2(1:2),'-o')
if length(CueType)>2
    hold on
    plot(CueType(3:4),pcorrCueT2(3:4),'-o')
end
ylim([0 1])
xlim([200 209])
xlabel('Cue (name)')
ylabel('Percent Correct')
title(b)
exportgraphics(gcf,['CorrectbyDir',b,'.png'],'Resolution',300)

corrTrials = data.ecodes.data(:,34);
actTarg = data.ecodes.data(:,35);
actTargPlot = actTarg==135&~isnan(actTarg);
figure
plot(movmean(corrTrials,25))
hold on
plot(actTargPlot,'s')
pbaspect([4 1 1])
xlabel('Trials')
ylabel('Percent Correct')
title(['Moving Avg. ',b])
% saveas(gcf,['MovAvg',b],'png')
exportgraphics(gcf,['MovAvg',b,'.png'],'Resolution',400)

sacEndX = data.ecodes.data(:,43);
sacEndY = data.ecodes.data(:,44);

figure 
subplot(2,1,1)
plot(sacEndX(corrTrials==1),sacEndY(corrTrials==1),'.b')
title(['Sacc End Correct ',b])
subplot(2,1,2)
plot(sacEndX(corrTrials==0),sacEndY(corrTrials==0),'.r')
xlabel('x eye position')
ylabel('y eye position')
title(['Sacc End Incorrect ',b])
saveas(gcf,['SaccEnd',b],'png')


figure
centerInd = find(data.ecodes.data(:,38)==180);
for t = 1:length(centerInd)
    subplot(2,1,corrTrials(centerInd(t))+1)
    eyeX = data.analog.data(centerInd(t),2).values(1200:(end-800));
    eyeY = data.analog.data(centerInd(t),3).values(1200:(end-800));
    hold on
    plot(eyeX,eyeY)
end

figure
notcenterInd = find(data.ecodes.data(:,38)~=180);
for t = 1:length(notcenterInd)
    subplot(2,1,corrTrials(notcenterInd(t))+1)
    eyeX = data.analog.data(notcenterInd(t),2).values(1200:(end-800));
    eyeY = data.analog.data(notcenterInd(t),3).values(1200:(end-800));
    hold on
    plot(eyeX,eyeY)
end