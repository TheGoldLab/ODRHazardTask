%Get the file
addpath(genpath('/Users/lab/Documents/Alice/ODRHazardTask/Analysis/TrackTraining'))
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = 17;

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
CueAngCon = CueAng;
if sum(CueAng>270)>0
    CueAngLeftConv = CueAng(CueAng>270)-360;
    CueAngCon(CueAng>270)= CueAngLeftConv;
end
figure
plot(CueAngCon,pcorrCue,'o')
ylim([0 1])
xlabel('Cue (deg)')
ylabel('Percent Correct')
title(['Session Correct ', b])
% saveas(gcf,['CorrectbyCue',b],'png')
exportgraphics(gcf,['CorrectbyCue',b,'.png'],'Resolution',300)

CueType = unique(data.ecodes.data(:,30));
for c = 1:length(CueType)
   cueInd = data.ecodes.data(:,30)==CueType(c)&data.ecodes.data(:,29)==2;
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
title(['T2 correct ',b])
exportgraphics(gcf,['CorrectbyDir',b,'.png'],'Resolution',300)

corrTrials = data.ecodes.data(:,34);
actTarg = data.ecodes.data(:,35);
actTargPlot = (actTarg==135|actTarg==45)&~isnan(actTarg);
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

centerInd = find(data.ecodes.data(:,38)==0|data.ecodes.data(:,38)==180);
corrTrials = data.ecodes.data(:,34);
corrTrialsCenter = corrTrials(centerInd);
actTarg = data.ecodes.data(:,35);
actTargPlot = (actTarg==135|actTarg==45)&~isnan(actTarg);
actTargPlotCenter=actTargPlot(centerInd);
figure
plot(movmean(corrTrialsCenter,25))
hold on
plot(actTargPlotCenter,'s')
pbaspect([4 1 1])
xlabel('Trials')
ylabel('Percent Correct')
title(['Moving Avg. Center Cue ',b])
% saveas(gcf,['MovAvg',b],'png')
exportgraphics(gcf,['MovAvgCenter',b,'.png'],'Resolution',400)



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
% for t = 1:length(centerInd)
%     subplot(2,1,corrTrials(centerInd(t))+1)
%     eyeX = data.analog.data(centerInd(t),2).values(1200:(end-800));
%     eyeY = data.analog.data(centerInd(t),3).values(1200:(end-800));
%     hold on
%     plot(eyeX,eyeY)
% end
% 
% figure
% notcenterInd = find(data.ecodes.data(:,38)~=180);
% for t = 1:length(notcenterInd)
%     subplot(2,1,corrTrials(notcenterInd(t))+1)
%     eyeX = data.analog.data(notcenterInd(t),2).values(1200:(end-800));
%     eyeY = data.analog.data(notcenterInd(t),3).values(1200:(end-800));
%     hold on
%     plot(eyeX,eyeY)
% end

finTrialInd = find(~isnan(data.ecodes.data(:,13)));
sampfreq=data.analog.acquire_rate(end);
tTargAcq = data.ecodes.data(finTrialInd,13)-data.ecodes.data(finTrialInd,5);
for tr = 1:length(finTrialInd)
    acqX(tr) = data.analog.data(finTrialInd(tr),2).values(round(tTargAcq(tr)));
    acqY(tr) = data.analog.data(finTrialInd(tr),3).values(round(tTargAcq(tr)));
end
figure
plot(acqX,acqY,'.')
xlabel('x eye position')
ylabel('y eye position')
title(['Target Acquisition ',b])
saveas(gcf,['TargAcq',b],'png')


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
figure
subplot(1,3,1)
plot(CueChoose(:,1),CueChoose(:,2),'bo-')
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
plot(CueChoose(:,1),CueChoose(:,2),'ro-')
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
plot(CueChoose(:,1),CueChoose(:,2),'go-')
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
plot(CueChoose(:,1),CueChoose(:,2),'co-')
title('No change')
exportgraphics(gcf,['pChooseUp',b,'.png'],'Resolution',400)