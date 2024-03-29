function [stats_, badTrials_] = plotFIRA_rasterPSTH(trials, Lsb, spikei, ...
    raster_begin, raster_end, rate_begin, rate_end, markers, ...
    bin_size, sortras, ras_axs, ras_hs, psth_axs, psth_hs)
%function stats_ = plotFIRA_rasterPSTH(trials, Lsb, spikei, ...
%    raster_begin, raster_end, rate_begin, rate_end, markers, ...
%    bin_size, sortrast, ras_axs, ras_hs, psth_axs, psth_hs);
%
% Returns:
%   stats_ ... array (one per plot) of
%               mean sem std median iqr n
global FIRA

% default
stats_ = [];

% easy check for now
if nargin < 11 || isempty(FIRA) || isempty(FIRA.spikes.data) || ...
        isempty(trials) || isempty(Lsb) || isempty(spikei) || ...
        spikei < 1 || spikei > size(FIRA.spikes.data,2)
    return
end

% length of ras_hs & psth_hs are either 0 or same
% as the number of columns of Lsb
num_axes = size(Lsb, 2);

% stats_ contains for each axis information about the rate:
% mean sem std median iqr n
stats_ = nans(num_axes, 6);

%% check times...
Lbad = false(size(trials));

% raster_begin determines start time of the plotted raster
if isempty(raster_begin)
    raster_begin = zeros(size(trials));
else
    Lbad = Lbad | isnan(raster_begin);
end

% raster_end determines end time of the plotted raster
if isempty(raster_end)
    raster_end = zeros(size(trials));
else
    Lbad = Lbad | isnan(raster_end);
end

% rate_begin determines rate computation, for polar plot
% also used as "wrt" for raster plot
if isempty(rate_begin) % use for wrt
    rate_begin = zeros(size(trials));    
    wrt        = zeros(size(trials));
else
    wrt        = rate_begin;
    Lbad       = Lbad & isnan(rate_begin);
end

% rate_end is end of rate computation
if isempty(rate_end)
    rate_end = zeros(size(trials));
else
    Lbad = Lbad | isnan(rate_end);
end

% trim data
if any(Lbad)
    Lgood        = ~Lbad;
    trials       = trials(Lgood);
    Lsb          = Lsb(Lgood, :);
    markers      = markers(Lgood, :);
    raster_begin = raster_begin(Lgood);
    raster_end   = raster_end(Lgood);
    rate_begin   = rate_begin(Lgood);
    wrt          = wrt(Lgood);
    
    if nargout == 2
        badTrials_ = sum(Lbad);
    end
else
    if nargout == 2
        badTrials_ = 0;
    end
end

% keep track of axis limits (height, left, right)
lims = zeros(num_axes, 4);

%%%
% DATA LOOP
%%%
%
% loop through the axes (angles), plotting rasters/psths and saving stats
for i = 1:num_axes

    % get the trial indices
    trs    = find(Lsb(:,i));    
    num_tr = length(trs);
    raster = [];
    rate   = zeros(num_tr, 1);

    % sort, if neccesary
    if sortras == 2 % sort by raster time
        [I,Y] = sort(raster_end(trs) - raster_begin(trs));
        trs = trs(Y);
        
    elseif sortras == 3 % sort by rate time
        [I,Y] = sort(rate_end(trs) - rate_begin(trs));
        trs = trs(Y);
    end
    
    % (hopefully) faster way of keeping track of the number of
    %   trials per bin, for the psth
    mint = floor(min(raster_begin(trs)-wrt(trs)));
    maxt = ceil(max(raster_end(trs)-wrt(trs)));
    lent = maxt - mint;
    if ~isempty(psth_axs) && ~isempty(lent) && lent < 10000 && ...
            all(~isnan(raster_begin(trs))) && all(~isnan(raster_end(trs))) && ...
            all(~isnan(wrt(trs)))
        psth_times = zeros(1, lent+1);
        psth_bts   = floor(raster_begin(trs) - wrt(trs)) - mint + 1;
        psth_ets   = ceil(raster_end(trs) - wrt(trs)) - mint + 1;
    else
        lent       = [];
        psth_times = [];
    end

    
    % loop through the trials
    for ti = 1:num_tr

        t = trs(ti);

        % get spikes from this trial
        sp = FIRA.spikes.data{trials(t), spikei};

        if ~isempty(sp)

            % raster
            msp = sp(sp>=raster_begin(t)&sp<=raster_end(t));
            if ~isempty(msp)
                raster = [raster; msp-wrt(t) ti*ones(size(msp))];
            end

            % rate
            if rate_end(t) > rate_begin(t)
                rate(ti) = 1000*sum(sp>=rate_begin(t)&sp<=rate_end(t)) / ...
                    (rate_end(t)-rate_begin(t));
            end

            % psth
            if ~isempty(psth_hs) && ~isempty(lent)
                
                % this is the way we like to compute it -- increment counts
                % for the given time interval
                
                psth_times(psth_bts(ti):psth_ets(ti)) = ...
                    psth_times(psth_bts(ti):psth_ets(ti)) + 1;
                
            elseif ~isinf(raster_begin(t)) && ~isinf(raster_end(t))
                
                % this is a much slower way -- keep track of all the times
                psth_times = [psth_times; [raster_begin(t):raster_end(t)]'-wrt(t)];

            else
                
                % this way sucks -- we don't really even know when the
                % trial begins & ends
                psth_times = [psth_times; [min(sp):max(sp)]'-wrt(t)];
                
            end
        end
    end

    if ~isempty(raster)

        % compute stats
        stats_(i, :) = [mean(rate) std(rate) std(rate) ...
            median(rate,1) iqr(rate) length(rate)];

        % plot raster, markers
        if ~isempty(ras_hs)
            set(ras_hs(i, 1), 'XData', raster(:,1), 'YData', raster(:,2));
            for m = 1:size(markers,2)
                set(ras_hs(i,2+m), 'XData', markers(trs,m)-wrt(trs), ...
                    'YData', 1:num_tr);
            end
        end

        % plot psth
        if ~isempty(psth_hs)
            if ~isempty(lent)
                xs = [mint:bin_size:maxt]';
                ys = histc(raster(:,1), xs)*1000; % ctl & jig removed bin_size scale term, 6/27/05
                ys = ys(:);
                yd  = zeros(bin_size, length(xs));
                len = min(length(psth_times), length(xs)*bin_size);
                yd(1:len) = psth_times(1:len);
                yd = sum(yd,1)';
            else
                xs = [min(psth_times):bin_size:max(psth_times)]';
                ys = histc(raster(:,1), xs)*1000/bin_size;
                ys = ys(:);
                yd = histc(psth_times, xs);
            end
            Lp      = yd>0;
            ys(Lp)  = ys(Lp)./yd(Lp);
            ys(~Lp) = 0;
            set(psth_hs(i,1), 'XData', xs, 'YData', ys);
            lims(i, 4) = max(ys);
        end

        % keep track of axis limits
        lims(i, 1:3) = [max(raster(:,2)) min(raster(:,1)) max(raster(:,1))];
    end
end

% clear extra markers
if ~isempty(ras_hs) && size(markers,2) < 3
    set(ras_hs(:, size(markers,2)+3:end), 'XData', [], 'YData', []);
end

if ~isempty(ras_axs) || ~isempty(psth_axs)
    minx = min(lims(:,2));
    maxx = max(lims(:,3));
    if minx == maxx
        maxx = minx + 1;
    end
end

% set the raster axes, draw line at 0
if ~isempty(ras_axs)

    % set axes
    nrows = max(lims(:,1));
    if nrows == 0
        nrows = 1;
    end
    set(ras_axs, 'XLim', [minx maxx], 'YLim', [0 nrows]);

    % draw line at zero
    set(ras_hs(:,2), 'XData', [0 0], 'YData', [0 nrows]);

    % set marker style for raster
    if nrows > 20
        set(ras_hs(:,1), 'Marker', '.', 'MarkerSize', 4);
    else
        set(ras_hs(:,1), 'Marker', '+', 'MarkerSize', 2);
    end
end

% set the raster axes, draw line at 0
if ~isempty(psth_axs)

    % set axes
    nrows = max(lims(:,4));
    if nrows == 0
        nrows = 1;
    end
    set(psth_axs, 'XLim', [minx maxx], 'YLim', [0 nrows]);

    % draw line at zero
    set(psth_hs(:,2), 'XData', [0 0], 'YData', [0 nrows]);
end
