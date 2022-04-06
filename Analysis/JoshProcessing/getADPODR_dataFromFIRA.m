function data_ = getADPODR_dataFromFIRA(fileName, monkey, behaviorOnly)
% function data_ = getADPODR_dataFromFIRA(fileName, monkey, behaviorOnly)
%

%% Parse args
%
% Get the file and put it in the global FIRA data struct so we can use 
%   the FIRA utilities
arguments
    fileName = 'MM_2022_03_21_7_00';
    monkey = 'MM';
    behaviorOnly = true;
end

% BehaviorOnly flag indicates that we don't get spike data and therefore
%   don't look in the "Sorted" directors for the raw data but rather
%   in the monkey-specific directory
if behaviorOnly
    if strcmp(monkey, 'MM')
        rawDataPath = 'C:/Users/alice/Box/GoldLab/Data/Physiology/AODR/MrM';
%         rawDataPath = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/MrM';
    else
        rawDataPath = 'C:/Users/alice/Box/GoldLab/Data/Physiology/AODR/Cicero';
    end
    convertedDataPath = 'C:/Users/alice/Box/GoldLab/Data/Physiology/AODR/BehaviorOnlyConverted';
    spikeList = [];
else
    % Otherwise get data from Sorted files
    rawDataPath = 'C:/Users/alice/Box/GoldLab/Data/Physiology/AODR/Sorted';
    convertedDataPath = 'C:/Users/alice/Box/GoldLab/Data/Physiology/AODR/SortedConverted';
    spikeList = 'all';
end

%% Check for file
%
% If fileName is an index, use it to get the name from the file list
if isempty(fileName)
    fileName = 1;
end
if isnumeric(fileName)    
    files = dir(fullfile(rawDataPath, '*.plx'));    
    fileName = files(min(max(length(files)-fileName+1,1),length(files))).name(1:end-4);
end

% Check if already converted
convertedFile = fullfile(convertedDataPath, fileName);
if isempty(dir([convertedFile '.mat']))
    bFile(fullfile(rawDataPath, [fileName '.plx']), [], ...
        'spmADPODR5', convertedFile, spikeList, 49:51, 0, 1, 0, []);
end

%% Open converted file
%
global FIRA
disp(fileName)
openFIRA(convertedFile);

%% Collect data from FIRA and put into a struct
%
data_ = struct('fileName', fileName, 'ecodes', [], 'timing', [], 'spikes', {{}});

% Use all non-broken fixation trials
score = getFIRA_ec([], 'score');
Lgood = score ~= -2;
nTrials = sum(Lgood);

%% Collect ecodes into a Table
data_.ecodes = table;

% These are values we can pull directly from FIRA
ecodeNames = {          ...
    'trial_num',        ... % within a session
    'task_id',          ... % 1=memory saccade; 2+=ADPODR
    'hazard',           ...
    'correct_target',   ... % Index of correct target (1 or 2)
    'sample_id',        ... % neg are close to T1, pos areclose to T2
    'score',            ... % 0=err, 1=correct, -3=sample, -1=ncerr
    'choice',           ... % 1=T1, 2=T2
    'RT',               ...
    'sac_endx',         ... 
    'sac_endy',         ... 
    't1_x',             ...
    't1_y',             ...
    't2_x',             ...
    't2_y',             ...
    'sample_x',         ...
    'sample_y',         ... 
    };
for nn = 1:length(ecodeNames)
    data_.ecodes.(ecodeNames{nn}) = getFIRA_ec(Lgood, ecodeNames{nn});
end

% These are values we need to compute
computeNames = {        ...    
    'session',          ... % To fill in elsewhere
    'tacp',             ... % Trials after change point
    'llr_for_switch',   ... % recode LLR wrt stay/switch
    'choice_switch',    ... % Re: **PREVIOUS TRUE TARGET**
	'previous_reward'   ... % 0=none/error, 1=yes/correct)
    'correct_shown',    ... % 0/1 if correct target was shown
    };
for nn = 1:length(computeNames)
    data_.ecodes.(computeNames{nn}) = zeros(sum(Lgood),1);
end

%  Trials after change point
noChangeTrials = [false; diff(data_.ecodes.correct_target)==0];
counter = 0;
for ii = 1:nTrials
    if noChangeTrials(ii)
        counter = counter + 1;
    else
        counter = 0;
    end
    data_.ecodes.tacp(ii) = counter;
end

% LLR for switch
% Get LLRs encoded as favoring T1 (-) vs T2 (+)
data_.ecodes.llr_for_switch = getFIRA_ec(Lgood, 'llr');
% Convert all to favoring current T (+) vs not (-)
Lct = data_.ecodes.correct_target==1;
data_.ecodes.llr_for_switch(Lct) = -data_.ecodes.llr_for_switch(Lct);
% Flip signs on non-cp trials
Lncp = data_.ecodes.tacp~=0;
data_.ecodes.llr_for_switch(Lncp) = -data_.ecodes.llr_for_switch(Lncp);

%  Choice switch (0/1)
prevState = [nan; data_.ecodes.correct_target(1:end-1)];
data_.ecodes.choice_switch(data_.ecodes.choice~=prevState)=1;
data_.ecodes.choice_switch(~ismember(data_.ecodes.choice, [1 2]) | ...
    ~ismember(prevState,[1 2])) = nan;

%  Previous reward (0=none/error, 1=yes/correct)
OLscore = getFIRA_ec(Lgood, 'online_score'); % get "on-line" scores
data_.ecodes.previous_reward = [0; OLscore(1:end-1)==1];

%  Correct Target Shown? 0=no, 1=yes
mean_on = getFIRA_ec(Lgood, 'mean_on');
target_off = getFIRA_ec(Lgood, 'target_off');
data_.ecodes.correct_shown = ~isnan(mean_on) & (mean_on < target_off);

%% Collect trial-specific timing values
% Keeping this separate from the ecode matris just to make them
%   both easier to parse/deal with
data_.timing = table;

% These are values we can pull directly from FIRA
timingNames = { ...
    'ring_on', 'sample_on', 'mean_on', 'target_off', 'fp_off', ...
    'choices_on', 'fdbk_on', 'fdbk_off', 'targ_acq', 'all_off'};
for nn = 1:length(timingNames)
    data_.timing.(timingNames{nn}) = getFIRA_ec(Lgood, timingNames{nn});
end

% These are values we need to compute
data_.timing.sac_on = getFIRA_ec(Lgood, 'fp_off') + getFIRA_ec(Lgood, 'RT');

%% Collect spikes
if isfield(FIRA, 'spikes')
    data_.spikes = FIRA.spikes.data(Lgood, :);
    data_.spikeid = FIRA.spikes.id;
end