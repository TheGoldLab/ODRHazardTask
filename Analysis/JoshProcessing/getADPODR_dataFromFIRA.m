function data_ = getADPODR_dataFromFIRA(fileName)
%
% Get data for behavioral analyses
global FIRA

% Get the file and put it in the global FIRA data struct so we can use 
%   the FIRA utilities
if nargin < 1 || isempty(fileName)    
    fileName = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/Converted/MM_2022_02_22_5_44-sort01-3units';
end
[~, fn] = fileparts(fileName);
openFIRA(fileName);

%% Collect data from FIRA and put into a struct:
data_ = struct('fileName', fn, 'ecodes', [], 'timing', [], 'spikes', {{}});

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
    'sac_endx',         ... % saccade end x
    'sac_endy',         ... % saccade end y
    't1_x',             ... % Target 1 x
    't1_y',             ... % Target 1 y
    't2_x',             ... % Target 2 x
    't2_y',             ... % Target 2 y
    'sample_x',         ... % Sample x
    'sample_y',         ... % Sample y
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

%  23. Correct Target Shown? 0=no, 1=yes
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
    'choices_on', 'fdbk_on', 'fdbk_off', 'all_off'};
for nn = 1:length(timingNames)
    data_.timing.(timingNames{nn}) = getFIRA_ec(Lgood, timingNames{nn});
end

% These are values we need to compute
data_.timing.sac_on = getFIRA_ec(Lgood, 'fp_off') + getFIRA_ec(Lgood, 'RT');

%% Collect spikes
data_.spikes = FIRA.spikes.data(Lgood, :);