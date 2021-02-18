function [flag_noTrial] = buildFIRA_cleanup(call_spm)
% function buildFIRA_cleanup(call_spm)
%
%   calls dataType cleanup methods
%
% Argument:
%   call_spm ... flag to force call to spm 'cleanup' routine

% created by jig 11/01/04

global FIRA

% YF 20160712
    if isempty(FIRA.raw.trial.start_time)
        flag_noTrial = 1;
        return
    else
        flag_noTrial = 0;
    end

% cleanup
for fn = fieldnames(FIRA.spm)'
    cleanup(FIRA.spm.(fn{:}));
end

if nargin > 0 && call_spm && ~isempty(FIRA.header.spmF)
    feval(FIRA.header.spmF, 'cleanup');
end
