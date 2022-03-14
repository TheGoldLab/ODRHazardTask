function plotADPODR_recentBehaviorSummary(monk, sessions, rawDataPath, convertedDataPath)
% function plotADPODR_recentBehaviorSummary(monk, numSessions)
% 

% Check arguments
arguments
    monk = 'MM';
    sessions = 5;
    rawDataPath = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/Sorted';
    convertedDataPath = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/Converted';
end

% Check for raw files
rawFiles = dir(fullfile(rawDataPath, [monk '*.plx']));

% sessions determines which sessions to plot -- either:
%   scalar: n most recent sessions
%   array: indices of sessions from file list
%   cell array: list of file names
if isscalar(sessions)
    sessions = 1:min(sessions, length(rawFiles));
    if isempty(sessions)
        disp('NO DATA FOUND')
        return
    end
elseif isnumeric(sessions)
    sessions = sessions(sessions<=length(rawFiles));
elseif iscell(sessions)
    names = {rawFiles.name};
    inds  = [];
    for ss = 1:length(sessions)
        inds = cat(2, inds, ...
            find(strncmp(sessions{ss}, names, length(sessions{ss}))));
    end
    sessions = inds;
end
numSessions = length(sessions);
    
% Open the figure
wid     = 16; % total width
hts     = 16/numSessions.*0.9;
cols    = num2cell(ones(1,numSessions).*3);
[axs,~] = getPLOT_axes([], wid, hts, cols, 2, 1.5, [], 'ADPODR', true);
set(axs,'Units','normalized');

% Loop through the sessions
for ii = 1:numSessions
    
    % Get and print the filename (minus suffix)
    fileName = rawFiles(end-sessions(ii)+1).name(1:end-4);    
    disp(fileName)
    
    % Check if already converted
    convertedFile = fullfile(convertedDataPath, fileName);
    if isempty(dir([convertedFile '.mat']))
        
        % need to convert    
        bFile(fullfile(rawDataPath, fileName), [], ...
            'spmADPODR5', convertedFile, 'all', 49:51, 0, 1, 0, []);
    end
    
    % Call plot function with appropriate axes
    plotADPODR_sessionBehavior(fileName, axs((ii-1)*3+(1:3)));    
end

    
    