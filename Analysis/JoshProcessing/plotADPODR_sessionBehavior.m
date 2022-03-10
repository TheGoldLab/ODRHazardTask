function plotADPODR_sessionBehavior(fileName)
% function plotADPODR_sessionBehavior(fileName)

%% Get data
% ecode names are:
%   trial_num, task_id, hazard, correct_target, sample_id, score, 
%   choice, RT, sac_endx, sac_endy, t1_x, t1_y, t2_x, t2_y, sample_x, 
%   sample_y, session, tacp, llr_for_switch, choice_switch, 
%   previous_reward, correct_shown
if nargin < 1
    fileName = [];
end
data = getADPODR_dataFromFIRA(fileName);

%% Collect relevant data
% PMF
% Independent variable is "signed cues" -- which is cue location
%   signed by corresponding "LLR for switch"
hazards    = nonanunique(data.ecodes.hazard);
numHazards = length(hazards);
sCues      = sign(data.ecodes.llr_for_switch).*abs(data.ecodes.sample_id);
cues       = nonanunique(sCues);
numCues    = length(cues);

% Dependent variable is choice switch re: previous TRUE STATE
%   (note that this is not the same as a choice switch, because 
%   the previous choice might have been in error).
Lswitch    = data.ecodes.choice_switch==1;
Lstay      = data.ecodes.choice_switch==0;

% collect data in pdat
%   dim 1 is LLR
%   dim 2 is pmf, cmf T1/T2 choices
%   dim 3 is hazard
%   dim 4 is for data, ns, sem
pdat  = nans(numCues, 3, numHazards, 3); 
Lc1   = data.ecodes.choice==1;
Lc2   = data.ecodes.choice==2;
Ltask = data.ecodes.task_id>=2;
Lgood = Ltask & data.ecodes.score>=0;
% For each cue location
for ll = 1:numCues
    Lcue = Lgood & sCues==cues(ll);
    
    % For each hazard
    for hh = 1:numHazards
        Lh = Lcue & data.ecodes.hazard == hazards(hh);
        
        % collect data
        pdat(ll,:,hh,1) = [ ...
            sum(Lh&Lswitch)./sum(Lh&(Lswitch|Lstay)).*100, ...
            nanmean(data.ecodes.RT(Lh&Lc1)), ...
            nanmean(data.ecodes.RT(Lh&Lc2))];
        
        % collect ns
        pdat(ll,:,hh,2) = [ ...
            sum(Lh&(Lswitch|Lstay)), sum(Lh&Lc1), sum(Lh&Lc2)];

        % collect sems
        pdat(ll,2:3,hh,3) = [ ...
            nanse(data.ecodes.RT(Lh&Lc1)), ...
            nanse(data.ecodes.RT(Lh&Lc2))];
    end
end

%% Plotz
co = {[4 94 167]./255, [194 0 77]./255};
    
%% 1. PMF
% First get logistic fit
%   input matrix is h1_bias, h2_bias, llr for switch, choice switch
ldat = cat(2, zeros(sum(Lgood), 2), sCues(Lgood), data.ecodes.choice_switch(Lgood));
ldat(data.ecodes.hazard(Lgood)==hazards(1), 1) = 1;
ldat(data.ecodes.hazard(Lgood)==hazards(2), 2) = 1;
ldat = ldat(isfinite(ldat(:,4)),:);
fits = logist_fit(ldat, 'lu1');
xax = (-4:0.01:4)';
lx  = cat(2, ones(size(xax)), xax);

% show it
% Get max n to scale points
maxn = max(reshape(pdat(:,1,:,2),[],1));
subplot(4,1,1); cla reset; hold on;
plot([-4 4], [50 50], 'k:');
plot([0 0], [0 100], 'k:');
% For each hazard
for hh = 1:numHazards
    % Plot raw data per LLR
    for ll = 1:numCues
        if pdat(ll,1,hh,2) > 0
            plot(cues(ll), pdat(ll,1,hh,1), '.', 'MarkerSize', pdat(ll,1,hh,2)/maxn*50, ...
                'Color', co{hh});
        end
    end
    plot(xax, logist_valLU1(fits([hh 3 4]), lx).*100, '-', 'Color', co{hh}, 'LineWidth', 2);
    
    % bias is just log((1-H)/H) -- see Glaze et al 2015
    h=text(-3.5, 60+(hh-1)*15, sprintf('fit H=%0.2f (true=%.2f)', ...
        1./(1+exp(-fits(hh))), hazards(hh)));
    set(h, 'Color', co{hh}, 'FontSize', 20)
end
axis([-4 4 0 100]);
xlabel('Cue location(-=stay, +=switch)')
ylabel('Prob(switch)');
title(strrep(data.fileName, '_', '-'))
set(gca,'FontSize', 14);

%% 2 & 3. CMF T1 & T2 choices
MIN_N = 5;
for cc = 1:2
    
    % Set up plot
    subplot(4,1,cc+1); cla reset; hold on;
    plot([0 0], [0 1000], 'k:');
    
    % For each hazard
    for hh = 1:numHazards
        
        % Plot error bars, means
        Lg = pdat(:,cc+1,hh,2)>=MIN_N;
        xs = cues(Lg);
        ys = pdat(Lg,cc+1,hh,1);
        es = pdat(Lg,cc+1,hh,3);
        plot([xs xs]', [ys-es ys+es]', '-', 'Color', co{hh});
        plot(xs, ys, 'o-', 'MarkerSize', 10, 'Color', co{hh}, ...
            'MarkerFaceColor', co{hh});
    end
    axis([-4 4 100 300]);
    set(gca,'FontSize', 14);
    ylabel('RT (ms)')   
    title(sprintf('RT, choice=%d', cc))
end
xlabel('Cue location(-=stay, +=switch)')

%% 4. Eye positions
% Use different symbols for different target configurations
symbols = {'o', 's', 'd', 'v', '^', '<', '>'};
% Separate by choices: T1, T2, chose_cue, NC
colors = {[153 216 201]./255, [44 162 95]./255, 'm', 'k'}; 
Lchc = [data.ecodes.choice==1 data.ecodes.choice==2, ...
    data.ecodes.score==-3, data.ecodes.score==-1];
subplot(4,1,4); cla reset; hold on;
plot([-50 50], [0 0], 'k:');
plot([0 0], [-50 50], 'k:');
% Loop through each Target configuration
Ts = table2array(data.ecodes(:, {'t1_x', 't1_y', 't2_x', 't2_y'}));
uTs = unique(Ts, 'rows');
uTs = uTs(isfinite(uTs(:,4)),:);
% For each set of target positions
for rr = 1:size(uTs,1)
    
    % Get trials
    Ltrials = all(Ts==uTs(rr,:), 2);

    % For each choice type
    for cc = 1:size(Lchc,2)

        % Plot eye positions
        plot(data.ecodes.sac_endx(Ltrials&Lchc(:,cc)), ...
            data.ecodes.sac_endy(Ltrials&Lchc(:,cc)), ...
            symbols{rr}, 'MarkerSize', 5, ...
            'Color', colors{cc}, 'MarkerFaceColor', colors{cc});

        % plot target
        if cc <= 2
            plot(uTs(rr,(cc-1)*2+1), uTs(rr,(cc-1)*2+2), ...
                symbols{rr}, 'MarkerSize', 20, 'color', 'k', ...
                'MarkerFaceColor', 0.5.*ones(1,3));
            plot(uTs(rr,(cc-1)*2+1), uTs(rr,(cc-1)*2+2), ...
                symbols{rr}, 'MarkerSize', 10, ...
                'MarkerFaceColor', colors{cc});
        end
    end
end
axis([-30 30 -30 30]);
set(gca, 'FontSize', 14);
title('Eye positions')
xlabel('Horizontal position (dva)');
ylabel('Vertical position (dva)');
