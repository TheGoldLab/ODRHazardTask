%Put neural data in trial structure with behavior data
function postRecordNeuralCheck(data_path,save_path)
%Get the file
% addpath(genpath('/Users/lab/Documents/Alice/ODRHazardTask/Analysis/TrackTraining'))
% data_path = '/Users/lab/Desktop/MM_training/MatFiles';
% figPath = '/Users/lab/Desktop/MM_training/Figures';
% data_path = '/Users/lab/Desktop/Ci_training/MatFiles';
% figPath = '/Users/lab/Desktop/Ci_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = nfiles;

theFile = files(currentFile).name;
[~,b,~]=fileparts(theFile);
b(strfind(b,'_'))='-';

temp=load(theFile);
data=temp.data;
clear temp


cd(save_path)

global FIRA
FIRA = data;
clear data
getF = @(t,n) FIRA.ecodes.data(t, strcmp(n, FIRA.ecodes.name));

unitsInd = FIRA.spikes.unit>0;
trialStartTS = getF(':', 'trial_begin');
trialEndTS = getF(':', 'trial_end');

trialNum = getF(':', 'trial_num');

corr = getF(':', 'Offline_Score')==1;

finTrialInd = ~isnan(getF(':', 'targacq'));
aODRind = getF(':', 'taskid') == 2|getF(':', 'taskid') == 3;
VMGSind = getF(':', 'taskid') == 1;

if sum(VMGSind)>0
    VMGStrialNum = trialNum(VMGSind);
    unitsIndNum = find(unitsInd);
    for u = 1:sum(unitsInd)
        figure
        hold on
        for t = 1:length(VMGStrialNum)
            plot(FIRA.spikes.data{VMGStrialNum(end+1-t),unitsIndNum(u)},VMGStrialNum(end+1-t).*ones(size(FIRA.spikes.data{VMGStrialNum(end+1-t),...
            unitsIndNum(u)})),'.k')
        
        end
        title(FIRA.spikes.id(unitsIndNum(u)))
    end
    
end

if sum(aODRind)>0
    aODRtrialNum = trialNum(aODRind);
    unitsIndNum = find(unitsInd);
    for u = 1:sum(unitsInd)
        figure
        hold on
        for t = 1:length(aODRtrialNum)
            plot(FIRA.spikes.data{aODRtrialNum(end+1-t),unitsIndNum(u)},aODRtrialNum(end+1-t).*ones(size(FIRA.spikes.data{aODRtrialNum(end+1-t),...
            unitsIndNum(u)})),'.k')
        
        end
        title(FIRA.spikes.id(unitsIndNum(u)))
    end
    
end


end
