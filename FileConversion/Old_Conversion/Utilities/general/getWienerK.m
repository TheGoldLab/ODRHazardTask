function wk_ = getWienerK(in_vec, out_vec, max_lag)
% function wk_ = getWienerK(in_vec, out_vec, max_lag)
%
% Written 3/18/04 by Joshua Gold

% check args
if nargin < 2 || isempty(in_vec) || isempty(out_vec)
    return;
end

if nargin < 3 || isempty(max_lag)
    max_lag = 10;
end

% get autocorrelation vector of input
ac = xcorr(in_vec, max_lag);

% get xcorr of input/output
xc = xcorr(in_vec, out_vec, max_lag);

% compute wiener kernel
wk_ = toeplitz(ac(max_lag+1:-1:1))\xc(max_lag+1:-1:1);
