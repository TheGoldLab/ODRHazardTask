%Get the file
addpath(genpath('/Users/lab/Documents/Alice/ODRHazardTask/Analysis/TrackTraining'))
data_path = 'C:\Users\alice\Desktop\Ci_training';%'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Cicero';
figPath = 'C:\Users\alice\Desktop\Ci_training';
% data_path = '/Users/lab/Desktop/MM_training/MatFiles';
% figPath = '/Users/lab/Desktop/MM_training/Figures';
% data_path = '/Users/lab/Desktop/Ci_training/MatFiles';
% figPath = '/Users/lab/Desktop/Ci_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = nfiles;

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

global FIRA
FIRA = data;
clear data
getF = @(t,n) FIRA.ecodes.data(t, strcmp(n, FIRA.ecodes.name));

% Make some plots

if sum(unique(getF(':', 'taskid')) == 1)
    VMGSind = getF(':', 'taskid') == 1;
    
    sacEndX = getF(VMGSind, 'sac_endx');
    sacEndY = getF(VMGSind, 'sac_endy');
    
    figure
    plot(sacEndX,sacEndY,'.')
    xlim([-30 30])
    ylim([-30 30])
    title(['VGS/MGS Saccades ',b])
    saveas(gcf,['VSGMGSSacc',b],'png')
end


if sum(unique(getF(':', 'taskid')) == 0)
    calibind = getF(':', 'taskid') == 0;
    
    figure
    for t = 1:length(calibind)

        eyeX = FIRA.analog.data(t,2).values;
        eyeY = FIRA.analog.data(t,3).values;
        hold on
        plot(eyeX,eyeY)
    end
    xlim([-30 30])
    ylim([-30 30])
    title(['Calibration ',b])
    saveas(gcf,['Calib',b],'png')
end


if sum(unique(getF(':', 'taskid')) == 2)|...
   sum(unique(getF(':', 'taskid')) == 3)   

aODRind = getF(':', 'taskid') == 2|getF(':', 'taskid') == 3;

corrInd = getF(aODRind, 'Offline_Score') == 1;
errInd = getF(aODRind, 'Offline_Score') == 0;
ncInd = getF(aODRind, 'Offline_Score') == -1;
ccInd = getF(aODRind, 'Offline_Score') == -3;

sacEndX = getF(aODRind, 'sac_endx');
sacEndY = getF(aODRind, 'sac_endy');

actTarg = getF(aODRind, 'active_target');
t1_x = nanmedian(getF(':', 't1_x'));
t1_y = nanmedian(getF(':', 't1_y'));
t2_x = nanmedian(getF(':', 't2_x'));
t2_y = nanmedian(getF(':', 't2_y'));
% t1_a = nanmedian(getF(':', 't1_angle'));
% t2_a = nanmedian(getF(':', 't2_angle'));


CueAng = unique(getF(aODRind, 'sample_angle'));
for c = 1:length(CueAng)
   cueInd = getF(aODRind, 'sample_angle')==CueAng(c);
   pcorrCue(c) =  sum(corrInd&cueInd)./(sum(corrInd&cueInd)+sum(errInd&cueInd));   
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


CueType = unique(getF(aODRind, 'trialid'));
for c = 1:length(CueType)
   cueInd = getF(aODRind, 'trialid')==CueType(c)&getF(aODRind, 'taskid')==2;
   pcorrCueT2(c) =  sum(corrInd&cueInd)./(sum(corrInd&cueInd)+sum(errInd&cueInd));   
end
figure
plot(CueType(CueType<209),pcorrCueT2(CueType<209),'-o')

    hold on
    plot(CueType(CueType>=209&CueType<300),pcorrCueT2(CueType>=209&CueType<300),'-o')

ylim([0 1])
xlim([200 220])
xlabel('Cue (name)')
ylabel('Percent Correct')
title(['T2 correct ',b])
exportgraphics(gcf,['CorrectbyDir',b,'.png'],'Resolution',300)

% % corrTrials = data.ecodes.data(:,34);
% % corrTrials = data.ecodes.data(:,40)==1;
% % actTarg = data.ecodes.data(:,35);
actTargPlot = (actTarg==135|actTarg==45)&~isnan(actTarg);
actTargPlot= actTargPlot(corrInd|errInd);
figure
plot(movmean(corrInd(corrInd|errInd),25))
hold on
plot(actTargPlot,'s')
pbaspect([4 1 1])
xlabel('Trials')
ylabel('Percent Correct')
title(['Moving Avg. ',b])
% saveas(gcf,['MovAvg',b],'png')
exportgraphics(gcf,['MovAvg',b,'.png'],'Resolution',400)
% 
% centerInd = find(data.ecodes.data(:,38)==0|data.ecodes.data(:,38)==180);
% corrTrials = data.ecodes.data(:,34);
% corrTrialsCenter = corrTrials(centerInd);
% actTarg = data.ecodes.data(:,35);
% actTargPlot = (actTarg==135|actTarg==45)&~isnan(actTarg);
% actTargPlotCenter=actTargPlot(centerInd);
% figure
% plot(movmean(corrTrialsCenter,25))
% hold on
% plot(actTargPlotCenter,'s')
% pbaspect([4 1 1])
% xlabel('Trials')
% ylabel('Percent Correct')
% title(['Moving Avg. Center Cue ',b])
% % saveas(gcf,['MovAvg',b],'png')
% exportgraphics(gcf,['MovAvgCenter',b,'.png'],'Resolution',400)
% 
% % 
% % 
% % sacEndX = data.ecodes.data(:,43);
% % sacEndY = data.ecodes.data(:,44);
% 
% figure 
% subplot(2,1,1)
% plot(sacEndX(corrTrials==1),sacEndY(corrTrials==1),'.b')
% title(['Sacc End Correct ',b])
% subplot(2,1,2)
% plot(sacEndX(corrTrials==0),sacEndY(corrTrials==0),'.r')
% xlabel('x eye position')
% ylabel('y eye position')
% title(['Sacc End Incorrect ',b])
% saveas(gcf,['SaccEnd',b],'png')
% 
% 
% 
% 
% % figure
% % for t = 1:length(centerInd)
% %     subplot(2,1,corrTrials(centerInd(t))+1)
% %     eyeX = data.analog.data(centerInd(t),2).values(1200:(end-800));
% %     eyeY = data.analog.data(centerInd(t),3).values(1200:(end-800));
% %     hold on
% %     plot(eyeX,eyeY)
% % end
% % 
% % figure
% % notcenterInd = find(data.ecodes.data(:,38)~=180);
% % for t = 1:length(notcenterInd)
% %     subplot(2,1,corrTrials(notcenterInd(t))+1)
% %     eyeX = data.analog.data(notcenterInd(t),2).values(1200:(end-800));
% %     eyeY = data.analog.data(notcenterInd(t),3).values(1200:(end-800));
% %     hold on
% %     plot(eyeX,eyeY)
% % end

% finTrialInd = find(~isnan(data.ecodes.data(:,13)));
finTrialInd = find(~isnan(getF(':', 'targacq')));
aODRfindTrialInd = intersect(finTrialInd,find(aODRind));
sampfreq=FIRA.analog.acquire_rate(end);
% tTargAcq = data.ecodes.data(finTrialInd,13)-data.ecodes.data(finTrialInd,5);
tTargAcq = getF(aODRfindTrialInd, 'targacq')-getF(aODRfindTrialInd, 'fp_on');
for tr = 1:length(aODRfindTrialInd)
%     if intersect(finTrialInd(tr),find(aODRind))
        acqX(tr) = FIRA.analog.data(aODRfindTrialInd(tr),2).values(round(tTargAcq(tr)));
        acqY(tr) = FIRA.analog.data(aODRfindTrialInd(tr),3).values(round(tTargAcq(tr)));
%     else
%         acqX(tr) = nan;
%         acqY(tr) = nan;
%     end
end
figure
plot(acqX,acqY,'.')
xlabel('x eye position')
ylabel('y eye position')
title(['Target Acquisition ',b])
saveas(gcf,['TargAcq',b],'png')

figure
for t = 1:length(aODRfindTrialInd)
    
    eyeX = FIRA.analog.data(aODRfindTrialInd(t),2).values(1200:round(tTargAcq(t)));
    eyeY = FIRA.analog.data(aODRfindTrialInd(t),3).values(1200:round(tTargAcq(t)));
    hold on
    plot(eyeX,eyeY)
end


choiceDir = acqY>0;

cueFin =getF(intersect(finTrialInd,find(aODRind)), 'sample_angle');
actTargFin = getF(intersect(finTrialInd,find(aODRind)), 'active_target');
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
   taskInd = getF(intersect(finTrialInd,find(aODRind)), 'taskid') == 3;
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
   taskInd = (getF(intersect(finTrialInd,find(aODRind)), 'taskid') ==2)&(~actTargFinUp)&~isnan(actTargFin);
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
   taskInd = (getF(intersect(finTrialInd,find(aODRind)), 'taskid') ==2)&(actTargFinUp);
    cueInd = intersect(find(cueFin==CueAng(c)),find(taskInd));
   pChooseUp(c) =  sum(choiceDir(cueInd))./length(choiceDir(cueInd));   
end
CueChoose = [CueAngCon pChooseUp'];
CueChoose = sortrows(CueChoose);
CueChoose(isnan(CueChoose(:,2)),:)=[];
plot(CueChoose(:,1),CueChoose(:,2),'co-')
title('Low Hazard')
exportgraphics(gcf,['pChooseUp',b,'.png'],'Resolution',400)

end