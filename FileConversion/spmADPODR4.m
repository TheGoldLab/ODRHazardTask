function spmADPODR4(func)
% function spmLDdots(func)
%
% File describing how to interpret data for FIRA
% Use input argument "func" to switch functions:
%   'init'      ... fills FIRA.spm
%   'trial'     ... parse current trial (in FIRA.trial)
%   'cleanup'   ... whatever you need to do at the end
%
% Returns:
%   nada, but can change FIRA
%

% created by jig 10/21/04
% modified from Sharath's spm724rg.m
% LD 2-21-2007
% LD 9-17-08 added returnsac, omitrew Omitrew only for data after 9/15/08

% FOR DOTS TASK, CHOICE 1 IS DEFINED AS UP/LEFT; CHOICE2 IS DEFINED AS
% DOWN/RIGHT. NOTE THAT FOR REWARD ASYMMETRY TASKS, DEFINITION OF REWCONT=0
% MIGHT NOT MEAN CHOICE1 IS PAIRED WITH SMALL REWARD.

% edited by Yunshu 5-17-2016
% change how rewsize is computed given more than one beep
% now the rewsize is a little shorter than (rewoff-rewon), because rew on is
% the first time the tti pulse for valve opening (1st in the four pulses) is sent, and rew off is
% the last time the tti pulse for valve closing (3rd in the four pulses) is sent, which means the
% difference includes the time of interpulse intervel.

% Yunshu 5-23-2016
% Before 20150202 ntti should be 2n, after that, ntti should be 4n. If
% ntti is other strange number, set rewsize to 0 and rew off = rew on.



global FIRA
declareEC_ecodes_ADPODR2;

%% if called with empty FIRA.spm, fill it and return
if strcmp(func, 'init')
    declareEC_ecodes_ADPODR2;
    
    % useful ecode markers for tagged ecodes
    cb = 7000;
    cm = 6000;
    cx = 8000;
    
    % FIRA.spm is a struct with fields corresponding to "dataTypes"
    % that have values that are argument lists to the given dataType
    % constructor method (empty means use defaults).
    % (or, if you prefer, an nx2 cell array, with the first column
    % corresponding to the fieldnames, the second to the arglist.)
    %
    % Thus, to create FIRA with dataType "ecodes" including a list of
    % computed values called "one" and "two", use:
    % spm.ecodes = {'compute', {'one' 'value'; 'two' 'value}}
    %
    % buildFIRA_init then converts these descriptors into actual dataType
    % objects used to parse data. Note that the arglists to buildFIRA_init
    % can override the choices of dataTypes here (so you can use the same
    % base spm file to create different FIRAs)
    %
    % see FIRA/build/dataTypes for possible dataTypes
    %     FIRA/build/buildFIRA_init.m for details of how the spm struct is
    %       used
    FIRA.spm = struct( ...
        'trial',  {{ ...
        'startcd', 1005,                             ...
        'anycd',   [EC.CORRECTCD EC.WRONGCD EC.NOCHCD EC.BRFIXCD], ...
        'allcd',   [] ...
        }}, ...
        'ecodes', {{  ...
        'extract', {   ...  % values to extract & save in FIRA{2}: <name> <type> <method> <args>
        
        % Codes about Timing of task:
        'fp_on'         'time'      @getEC_time     {EC.FPONCD      1}; ...  %FP on
        'ring_on'       'time'      @getEC_time     {EC.RINGON      1};...
        
        'sample_on'     'time'      @getEC_time     {EC.TRGC1CD    1}; ...   %sample on
        'mean_on'       'time'      @getEC_time     {EC.TARGONCD 1}; ...
        'graphics_off (Mem)'   'time' @getEC_time	{EC.TARGOFFCD   1}; ...
        'fp_off'        'time'      @getEC_time     {EC.FPOFFCD     1}; ...
        'choices_on'	'time'      @getEC_time     {EC.TRGC3CD     1};...
        'sac_on'        'time'      @getEC_time     {EC.SACMAD      1}; ...
        'targacq'       'time'      @getEC_time     {EC.TRGACQUIRECD    1}; ...
        'fdbkon'        'time'      @getEC_time     {EC.FDBKONCD    1}; ...
        'fdbkoff'       'time'      @getEC_time     {EC.FDBKOFFCD    1}; ...
        'All_Off'       'time'      @getEC_time     {EC.ALLOFFCD   1}; ...
        
        %Codes about correct or not, both correct target and sample are
        %considered 'correct' here'
        'OLcorrect'     'time'      @getEC_time     {EC.CORRECTCD    1}; ...
        'OLerror'       'time'      @getEC_time     {EC.WRONGCD    1}; ...
        'OLncerr'       'time'      @getEC_time     {EC.NOCHCD    1}; ...
        'OLbrfix'       'time'      @getEC_time     {EC.BRFIXCD    1}; ...
        %Location, id, and environment info
        'fp_x'          'value'     @getEC_tagged	{EC.I_FIXXCD    cb    cm    cx  0.1};   ...
        'fp_y'          'value'     @getEC_tagged	{EC.I_FIXYCD    cb    cm    cx  0.1};   ...
        't1_x'          'value'     @getEC_tagged	{EC.I_TRG1XCD   cb    cm    cx  0.1};   ...
        't1_y'          'value'     @getEC_tagged   {EC.I_TRG1YCD   cb    cm    cx  0.1};   ...
        'sample_x'      'value'     @getEC_tagged   {EC.I_TRG3XCD   cb    cm    cx  0.1};   ...
        'sample_y'      'value'     @getEC_tagged   {EC.I_TRG3YCD   cb    cm    cx  0.1};   ...
        't2_x'          'value'     @getEC_tagged   {EC.I_TRG2XCD   cb    cm    cx  0.1};   ...
        't2_y'          'value'     @getEC_tagged   {EC.I_TRG2YCD   cb    cm    cx  0.1};   ...
        'taskid'        'id'        @getEC_tagged   {EC.I_TASKIDCD  cb    cm    cx  1};     ...
        'trialid'       'id'        @getEC_tagged   {EC.I_TRIALIDCD cb    cm    cx  1};     ...
        'seed_base'     'value'     @getEC_tagged   {EC.I_DTVARCD   cb    cm    cx  1};     ...
        'Haz'           'value'     @getEC_tagged   {EC.I_T1IDH   cb    cm    cx  .001};     ...
        'Noise'     	'value'     @getEC_tagged   {EC.I_T1SIGMA   cb    cm    cx  .1};     ...
        
        }, ...
        'compute', { ...                % names/types of fields to compute & save later in FIRA{2}
        'rew_on'        'id';       ...
        'active_target'	'id';       ...
        't1_angle'      'id';       ...
        't2_angle'      'id';       ...
        'sample_angle'	'id';       ...
        'sac_angle'     'id';       ...
        'Offline_Score' 'id';       ... % offline score: 1=correct, 0=error, -1=nc, -2=brfix,-3=sample
        'OLscore'       'value';    ... % online score: 1=correct, 0=error, -1=nc, -2=brfix
        'scorematch'	'value';    ...
        'sac_endx'      'value';    ...
        'sac_endy'      'value';    ...
        'sac_lat'       'value';    ...
        'sac_dur'       'value';    ...
        'sac_vmax'      'value';    ...
        'sac_amp'       'value';    ...
        'sac_on_offline' 'value';   ...
        'choice'        'id';       ...
        }, ...
        'tmp', { ...  % values to extract & use temporarily: <name> <type> <method> <args>
        }}}, ...
        'spikes',  [], ...
        'analog',  [], ...
        'dio',     []);
    
    return
    
    %% parse the trial ... check whether anything is given in FIRA.trial
elseif strcmp(func, 'trial')
    
    % note:
    %   for recordings before 4/22/09, All_OffCD is the feedback for both correct and error trials.
    %     Reward follows All_OffCD in correct trials
    %   for recordings afterward,
    %     in error trials, FDBKONCD is dropped 400 ms after TRGACQUIRECD, All_OffCD dropped 400 ms later.
    %     in correct trials, REW DIO dropped 800 ms after TRGACQUIRECD
    %             if (getFIRA_ec(tti, 'taskid')==3)
    %                 setFIRA_ec(tti, {'fdbkon'}, getFIRA_ec(tti, {'All_Off'}));
    %             end
    % get this trial index
    tti = FIRA.raw.trial.good_count;
    taskid = getFIRA_ec(tti, 'taskid');
    
    %     trialid = getFIRA_ec(tti, {'trialid'});
    % fixation off variable name
    if taskid>=1
        sfixoff='fp_off';
    else
        sfixoff='';
    end
    
    % Get reward time from dio
%     if isfield(FIRA, 'dio')
%         
%         %                 rewon = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
%         %                 rewoff = getDIO_time(FIRA.dio{tti}, 0, 2, 2);
%         %                 rewsize = rewoff - rewon;
%         
%         %                 Yunshu 2016-05-17
%         %                 deal with >1 beep rew by sum up all the valve opening time
%         fname = FIRA.header.filename{1,1};
%         fname_date = str2double(fname(5:10));
%         if fname_date > 150203
%             ndio_standard = 4;
%         else
%             ndio_standard = 2;
%         end
%         dio_tti = FIRA.dio{tti};
%         ntti = sum(dio_tti(:,2)==0 & dio_tti(:,3)==2);
%         n_beep = floor(ntti/ndio_standard);
%         if mod(ntti,ndio_standard)~=0
%             % sprintf('trial %.d mod(ntti)/%.d ~= 0',tti,ndio_standard)
%             if mod(ntti,ndio_standard)==1
%                 t_open = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
%                 t_close = t_open;
%             elseif mod(ntti,ndio_standard)>1
%                 t_open = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
%                 t_close = getDIO_time(FIRA.dio{tti}, 0, 2, 2);
%             end
%         else
%             t_open = nan(n_beep,1);
%             t_close = nan(n_beep,1);
%             for i_beep = 1:n_beep
%                 ind_open = i_beep*4-3;
%                 ind_close = i_beep*4-1;
%                 t_open(i_beep) = getDIO_time(FIRA.dio{tti}, 0, 2, ind_open);
%                 t_close(i_beep) = getDIO_time(FIRA.dio{tti}, 0, 2, ind_close);
%             end
%         end
%         % rewsize_beep = t_close - t_open;
%         % rewsize = sum(rewsize_beep);
%         
%         % for if numel(t_open) = 0 because of plx error or truncation
%         if numel(t_open) == 0
%             % OLscore = -3; % this trial will be excluded in later analyses
%             setFIRA_ec(tti, 'rew_on', 0);
%         else
%             setFIRA_ec(tti, 'rew_on', 1);
%             %                       rewon = t_open(1);
%             %                       rewoff = t_close(end); % rewsize < (rewoff- rewon), because take away the time between pulses
%             %                       setFIRA_ec(tti, {'rew_off'}, rewoff);
%             %                       setFIRA_ec(tti, {'reward'}, rewsize);
%         end
%     end
    
    % parse saccades
    if ~isempty(sfixoff)
        
        % inline functions to compute angles, distances in deg from xy
        ang_deg  = @(x,y) mod(atan2(y,x)*180/pi+360,360);        
        ang_diff = @(x,y) abs(rad2deg(angdiff(deg2rad(x),deg2rad(y))));
        % dst_deg = @(x,y) sqrt(x^2 + y^2);
        
        % Set t1/t2/sample angles
        setFIRA_ec(tti, 't1_angle', ang_deg(getFIRA_ec(tti, 't1_x'), getFIRA_ec(tti, 't1_y')));
        setFIRA_ec(tti, 't2_angle', ang_deg(getFIRA_ec(tti, 't2_x'), getFIRA_ec(tti, 't2_y')));
        if taskid > 1
            setFIRA_ec(tti, 'sample_angle', ang_deg(getFIRA_ec(tti, 'sample_x'), getFIRA_ec(tti, 'sample_y')));
        else
            setFIRA_ec(tti, 'sample_angle', getFIRA_ec(tti, 't1_angle'));
        end
        
        % Get and parse saccades
        [sacs, bf] = getFIRA_saccades(tti, 2, true, 'horiz_eye', 'vert_eye', sfixoff, 'fp_x', 'fp_y', 2000);
        if bf
            % fixbreak or other errors
            setFIRA_ec(tti, 'Offline_Score', -2);
        elseif isempty(sacs) || ~isfinite(sacs(1,1))
            % ncerr
            setFIRA_ec(tti, 'Offline_Score', -1);
        else
            
            % Found saccade, parse and store stats
            sacs = sacs(isfinite(sacs(:,1)),:);
            nSacs = size(sacs, 1);
            
            % Check for very quick second saccade, use if less<threshold
            while nSacs > 1 && sacs(2,1) < (sacs(1,1)+sacs(1,2)+30)
                sacs(2,1) = sacs(1,1);
                sacs(2,2) = sacs(2,1)-sacs(1,1)+sacs(2,2);
                sacs(1,:) = [];
                nSacs = nSacs - 1;
            end
            
            % Store saccade stats
            setFIRA_ec(tti, 'sac_lat',  sacs(1,1));
            setFIRA_ec(tti, 'sac_dur',  sacs(1,2));
            setFIRA_ec(tti, 'sac_vmax', sacs(1,3));
            setFIRA_ec(tti, 'sac_endx', sacs(1,5));
            setFIRA_ec(tti, 'sac_endy', sacs(1,6));
            setFIRA_ec(tti, 'sac_angle', ang_deg(sacs(1,5), sacs(1,6)));
            setFIRA_ec(tti, 'sac_amp',  sacs(1,7));
            setFIRA_ec(tti, 'sac_on_offline', getFIRA_ec(tti, sfixoff) + sacs(1,1));
            
            % For memory & visual & ring task, parse outcome
            if any(taskid == 2:5)

                % Pick the active target to see if correct
                if getFIRA_ec(tti, 'trialid')-200<9
                    correctTargetAngle = getFIRA_ec(tti, 't1_angle');
                    errorTargetAngle = getFIRA_ec(tti, 't2_angle');
                else
                    correctTargetAngle = getFIRA_ec(tti, 't2_angle');
                    errorTargetAngle = getFIRA_ec(tti, 't1_angle');
                end
                setFIRA_ec(tti, 'active_target', correctTargetAngle);
                
                % Check for good saccade based on length and angle
                setFIRA_ec(tti, 'Offline_Score', -1); % Default -- no choice
                amp = sqrt(getFIRA_ec(tti, 'sac_endx')^2 + getFIRA_ec(tti, 'sac_endy')^2 );
                ang = getFIRA_ec(tti, 'sac_angle');
                % if amp > 2.5 && amp < 12
                MIN_SACCADE_AMPLITUDE = 5;
                MAX_SACCADE_AMPLITUDE = 18;
                MAX_ANGULAR_DISTANCE  = 20;
                if amp > MIN_SACCADE_AMPLITUDE && amp < MAX_SACCADE_AMPLITUDE
                    if ang_diff(ang, correctTargetAngle) < MAX_ANGULAR_DISTANCE
                        setFIRA_ec(tti, 'Offline_Score', 1); % Correct!
                    elseif ang_diff(ang, errorTargetAngle) < MAX_ANGULAR_DISTANCE
                        setFIRA_ec(tti, 'Offline_Score', 0); % Error!
                    elseif ang_diff(ang, getFIRA_ec(tti, 'sample_angle')) < MAX_ANGULAR_DISTANCE
                        setFIRA_ec(tti, 'Offline_Score', -3); % Cue!
                    end
                end
                
                % PLOTZ -- set to true to see trial-by-trial eye
                % positions
                if false
                    
                    if n<9
                        Targ0X=getFIRA_ec(tti, {'t2_x'});
                        Targ0Y=getFIRA_ec(tti, {'t2_y'});
                    else
                        Targ0X=getFIRA_ec(tti, {'t1_x'});
                        Targ0Y=getFIRA_ec(tti, {'t1_y'});
                    end
                    co = {'b' 'm' 'k'};
                    xs = FIRA.analog.data(tti, 2).values;
                    ys = FIRA.analog.data(tti, 3).values;
                    ts = (0:length(xs)-1) + FIRA.analog.data(tti, 2).start_time - getFIRA_ec(tti, sfixoff);
                    subplot(2,1,1); cla reset; hold on;
                    plot([-2500 1500], [0 0], 'k:');
                    plot([0 0], [-10 10], 'k-');
                    plot(ts, xs, 'r-');
                    plot(ts, ys, 'g-');
                    for ss = 1:size(sacs, 1)
                        plot(sacs(ss,[1 1]), [-10 10], '-', 'Color', co{ss}, 'LineWidth', 2);
                    end
                    axis([-2500 1500 -8 8]);
                    
                    subplot(2,1,2); cla reset; hold on;
                    plot(xs(ts>0), ys(ts>0), '-', 'Color', [0.5 0.5 0.5]);
                    plot(TargX, TargY, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
                    plot(Targ0X, Targ0Y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
                    for ss = 1:size(sacs, 1)
                        plot(sacs(ss, 5), sacs(ss, 6), 'd', 'MarkerSize', 15, 'Color', co{ss}, 'MarkerFaceColor', co{ss});
                    end
                    plot([-10 10], [0 0], 'k:');
                    plot([0 0], [-10 10], 'k:');
                    axis([-10 10 -10 10]);
                    r = input('next')
                end
            end
        end
        
        % Compare to on-line score
        Lscore = isfinite(getFIRA_ec(tti, {'OLbrfix', 'OLncerr', 'OLerror', 'OLcorrect'}));
        if any(Lscore)
            setFIRA_ec(tti, 'OLscore', find(Lscore,1)-3); % convert to -2 -> 1
            setFIRA_ec(tti, 'scorematch', getFIRA_ec(tti, 'OLscore') == getFIRA_ec(tti, 'Offline_Score'));
        end
    end
end
