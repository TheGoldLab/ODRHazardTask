%Get the file
data_path = '/Users/lab/Desktop/MM_training';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = 4;

theFile= files(currentFile).name;
temp=load(theFile);
data=temp.data;

%make a save path to  new folder
[~,b,~]=fileparts(theFile);
mkdir(b)
savePath = b;
b(strfind(b,'_'))='-';

cd(savePath)

% Make some plots

CueType = unique(data.ecodes.data(:,38));
for c = 1:length(CueType)
   cueInd = data.ecodes.data(:,38)==CueType(c);
   pcorrCue(c) =  sum(data.ecodes.data(cueInd,34))./length(data.ecodes.data(cueInd,34));   
end
figure
plot(CueType,pcorrCue,'-o')
ylim([0 1])
xlabel('Cue (deg)')
ylabel('Percent Correct')
title(b)
% saveas(gcf,['CorrectbyCue',b],'png')
exportgraphics(gcf,['CorrectbyCue',b,'.png'],'Resolution',300)

corrTrials = data.ecodes.data(:,34);
actTarg = data.ecodes.data(:,35);
actTarg
actTargPlot = actTarg==max(actTarg);
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


% figure
% hold on
% centerInd = find(data.ecodes.data(:,38)==CueType(3));
% for t = 1:length(centerInd)
%     eyeX = data.analog.data(centerInd(t),2).values;
%     eyeY = data.analog.data(centerInd(t),3).values;
%     plot(eyeX,eyeY)
% end