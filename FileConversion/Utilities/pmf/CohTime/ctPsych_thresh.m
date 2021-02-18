function thresh_ = ctPsych_thresh(fun, fits, pctc, vtime)
% function thresh_ = ctPsych_thresh(fun, fits, pctc, vtime)
%
% computes threshold from general ctPsych function:
%   coherence needed for "pctc" % correct at "vtime"
%   viewing time. no lapse, no bias
% 
if nargin < 3 || isempty(pctc)
    pctc = 0.8161;
end

if nargin < 4 || isempty(vtime)
    vtime = 1;
end

% yeegads. Find initial values of fit == 0 (not given
%   data), assume these correspond to lapse & bias
%   terms, set to 0
in0       = feval(fun);
L0        = in0(:,1) == 0;
tfits     = fits;
tfits(L0) = 0;

GRAIN     = 0.001;
cohs      = (GRAIN:GRAIN:1)';
tdat      = [cohs vtime.*ones(size(cohs,1),1) ones(size(cohs,1))];
ps        = feval(fun, tfits, tdat);
ti        = find(ps>=pctc,1);
if isempty(ti)
    thresh_ = 1.0;
else
    thresh_ = cohs(ti);
end
