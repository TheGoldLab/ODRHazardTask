function spmADPODR3(func)
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
        'fp_on'     'time'      @getEC_time         {EC.FPONCD      1}; ...  %FP on
        'ring_on'   'time'      @getEC_time         {EC.RINGON      1};...
        
        'sample_on'    'time'      @getEC_time         {EC.TRGC1CD    1}; ...   %sample on 
        'mean_on'   'time'          @getEC_time         {EC.TARGONCD 1}; ...
        'graphics_off (Mem)'   'time'      @getEC_time         {EC.TARGOFFCD   1}; ...
        'fp_off'    'time'      @getEC_time         {EC.FPOFFCD     1}; ...
        'choices_on' 'time'     @getEC_time         {EC.TRGC3CD     1};...
         'sac_on'    'time'      @getEC_time         {EC.SACMAD      1}; ...
        'targacq'    'time'     @getEC_time         {EC.TRGACQUIRECD    1}; ...
        'fdbkon'    'time'      @getEC_time         {EC.FDBKONCD    1}; ...
        'fdbkoff'    'time'      @getEC_time         {EC.FDBKOFFCD    1}; ...
        'All_Off'    'time'     @getEC_time         {EC.ALLOFFCD   1}; ...

    %Codes about correct or not, both correct target and sample are
    %considered 'correct' here'
        'OLcorrect' 'time'      @getEC_time         {EC.CORRECTCD    1}; ...
        'OLerror'   'time'      @getEC_time         {EC.WRONGCD    1}; ...
        'OLncerr'   'time'      @getEC_time         {EC.NOCHCD    1}; ...
        'OLbrfix'   'time'      @getEC_time         {EC.BRFIXCD    1}; ...
    %Location, id, and environment info
        'fp_x'      'value'     @getEC_tagged       {EC.I_FIXXCD    cb    cm    cx  0.1};   ...
        'fp_y'      'value'     @getEC_tagged       {EC.I_FIXYCD    cb    cm    cx  0.1};   ...
        't1_x'      'value'     @getEC_tagged       {EC.I_TRG1XCD   cb    cm    cx  0.1};   ...
        't1_y'      'value'     @getEC_tagged       {EC.I_TRG1YCD   cb    cm    cx  0.1};   ...
        'sample_x'      'value'     @getEC_tagged       {EC.I_TRG3XCD   cb    cm    cx  0.1};   ...
        'sample_y'      'value'     @getEC_tagged       {EC.I_TRG3YCD   cb    cm    cx  0.1};   ...
        't2_x'      'value'     @getEC_tagged       {EC.I_TRG2XCD   cb    cm    cx  0.1};   ...
        't2_y'      'value'     @getEC_tagged       {EC.I_TRG2YCD   cb    cm    cx  0.1};   ...
        'taskid'    'id'        @getEC_tagged       {EC.I_TASKIDCD  cb    cm    cx  1};     ...
        'trialid'   'id'        @getEC_tagged       {EC.I_TRIALIDCD cb    cm    cx  1};     ...
        'seed_base' 'value'     @getEC_tagged       {EC.I_DTVARCD   cb    cm    cx  1};     ...
         'Haz'  'value'     @getEC_tagged       {EC.I_T1IDH   cb    cm    cx  .001};     ...
        'Noise'  'value'     @getEC_tagged       {EC.I_T1SIGMA   cb    cm    cx  .1};     ...

        }, ...
        'compute', { ...                % names/types of fields to compute & save later in FIRA{2}
        'rew_on'        'id';...
        'active_target'    'id';...
         't1_angle'       'id' ; ...     
        't2_angle'       'id';...
        'sample_angle',  'id';...
        'sac_angle'      'id';...
        'Offline_Score'       'id'; ... % offline score: 1=correct, 0=error, -1=nc, -2=brfix,-3=sample
        'OLscore'       'value'; ...    % online score: 1=correct, 0=error, -1=nc, -2=brfix
        'scorematch'    'value' ; ...
        

        'sac_endx'      'value' ; ...
        'sac_endy'      'value' ; ...
        'sac_lat'       'value' ; ...
        'sac_dur'       'value' ; ...
        'sac_vmax'      'value' ; ...
        'sac_amp'       'value' ; ...
        'sac_on_offline'       'value' ; ...
        }, ...
        'tmp', { ...  % values to extract & use temporarily: <name> <type> <method> <args>
        }}}, ...
        'spikes',  [], ...
        'analog',  [], ...
        'dio',     []);

    return

    % parse the trial ... check whether anything is given in FIRA.trial
elseif strcmp(func, 'trial')

    % get this trial index
    tti = FIRA.raw.trial.good_count;
    taskid = getFIRA_ec(tti, {'taskid'});
  
%     trialid = getFIRA_ec(tti, {'trialid'});
    % compute the saccade-related values as defined in init and determine
    % the rest

    if taskid>=1
        sfixoff='fp_off';
    else
        sfixoff='';
    end 


            
            if isfield(FIRA, 'dio')
                  
%                 rewon = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
%                 rewoff = getDIO_time(FIRA.dio{tti}, 0, 2, 2);
%                 rewsize = rewoff - rewon;

%                 Yunshu 2016-05-17
%                 deal with >1 beep rew by sum up all the valve opening time
                  fname = FIRA.header.filename{1,1};
                  fname_date = str2double(fname(5:10));
                  if fname_date > 150203
                      ndio_standard = 4;
                  else
                      ndio_standard = 2;
                  end
                  dio_tti = FIRA.dio{tti};
                  ntti = sum(dio_tti(:,2)==0 & dio_tti(:,3)==2);
                  n_beep = floor(ntti/ndio_standard);
                  if mod(ntti,ndio_standard)~=0
                     % sprintf('trial %.d mod(ntti)/%.d ~= 0',tti,ndio_standard)
                      if mod(ntti,ndio_standard)==1
                          t_open = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
                          t_close = t_open;
                      elseif mod(ntti,ndio_standard)>1
                          t_open = getDIO_time(FIRA.dio{tti}, 0, 2, 1);
                          t_close = getDIO_time(FIRA.dio{tti}, 0, 2, 2);
                      end
                  else
                      t_open = nan(n_beep,1);
                      t_close = nan(n_beep,1);
                      for i_beep = 1:n_beep
                          ind_open = i_beep*4-3;
                          ind_close = i_beep*4-1;
                          t_open(i_beep) = getDIO_time(FIRA.dio{tti}, 0, 2, ind_open);
                          t_close(i_beep) = getDIO_time(FIRA.dio{tti}, 0, 2, ind_close);
                      end
                  end
                  rewsize_beep = t_close - t_open;
                  rewsize = sum(rewsize_beep);
                  
                  % for if numel(t_open) = 0 because of plx error or
                  % truncation
                  if numel(t_open) == 0
                      OLscore = -3; % this trial will be excluded in later analyses
                    setFIRA_ec(tti, {'rew_on'}, 0);

                  else
%                       rewon = t_open(1);
%                       rewoff = t_close(end); % rewsize < (rewoff- rewon), because take away the time between pulses
                       setFIRA_ec(tti, {'rew_on'}, 1);
%                       setFIRA_ec(tti, {'rew_off'}, rewoff);
%                       setFIRA_ec(tti, {'reward'}, rewsize);
                  end
            end
    
    sacs_ = []; bf_ = [];
    if ~isempty(sfixoff)
%         [sacs_, bf_] = getFIRA_saccadesLD(tti, 2, false , 'horiz_eye',  'vert_eye', sfixoff, 'fp_x', 'fp_y', 5000); 
        [sacs_, bf_] = getFIRA_saccadesLD_SACMAD(tti, 2, true , 'horiz_eye',  'vert_eye', 'sac_on', 'fp_x', 'fp_y', 5000); 
        d_saconset = getFIRA_ec(tti, {'sac_on'}) - getFIRA_ec(tti, {sfixoff}) - 200;
        if ~isempty(sacs_)
            sacs_(:,1) = sacs_(:,1) + d_saconset;
        end
%         fpwin = 3;
%         [sacs_, bf_] = getFIRA_saccadesLD_fpwin(tti, 2, false , 'horiz_eye',  'vert_eye', sfixoff, 'fp_x', 'fp_y', 5000, fpwin); 
%         [sacs_, bf_] = getFIRA_saccades(tti, 2, false , 'horiz_eye', 'vert_eye', sfixoff, 'fp_x', 'fp_y', 5000); 
        t_eyeshift = 0; d_eyeshift = 0;
        t1_x = getFIRA_ec(tti, {'t1_x'});
        t1_y = getFIRA_ec(tti, {'t1_y'});
        t1_angle = mod(atan2(t1_y, t1_x)*180/pi+ 360,360);
       t2_x = getFIRA_ec(tti, {'t2_x'});
       t2_y = getFIRA_ec(tti, {'t2_y'});
       t2_angle = mod(atan2(t2_y, t2_x)*180/pi+ 360,360);   
           taskid = getFIRA_ec(tti, {'taskid'});
if taskid>1
       sample_x = getFIRA_ec(tti, {'sample_x'});
        sample_y = getFIRA_ec(tti, {'sample_y'});
       sample_angle = mod(atan2(sample_y, sample_x)*180/pi+ 360,360);
elseif taskid==1
       sample_angle = mod(atan2(t1_y, t1_x)*180/pi+ 360,360);

end
       
       
        
       setFIRA_ec(tti, {'t2_angle'}, t2_angle);
        setFIRA_ec(tti, {'t1_angle'}, t1_angle);
        setFIRA_ec(tti, {'sample_angle'}, sample_angle);
        
        if bf_
              %fixbreak or other errors
            setFIRA_ec(tti, {'Offline_Score'}, -2);
        elseif isempty(sacs_)
               %ncerr
            setFIRA_ec(tti, {'Offline_Score'}, -1);
        else
            [nSacs,m] = size(sacs_);
            
            
                      %  setFIRA_ec(tti, {'fix_time'}, getFIRA_ec(tti, {'eyefix'}) - getFIRA_ec(tti, {'fp_on'}));
            setFIRA_ec(tti, {'sac_lat'}, sacs_(1,1));
            setFIRA_ec(tti, {'sac_dur'}, sacs_(1,2));
            setFIRA_ec(tti, {'sac_vmax'}, sacs_(1,3));
            setFIRA_ec(tti, {'sac_endx'}, sacs_(1,5));
            setFIRA_ec(tti, {'sac_endy'}, sacs_(1,6));
        sac_x = getFIRA_ec(tti, {'sac_endx'});
        sac_y = getFIRA_ec(tti, {'sac_endy'});
        sac_angle = mod(atan2(sac_y, sac_x)*180/pi+ 360,360);
        
        setFIRA_ec(tti, {'sac_angle'}, sac_angle);
            
            
            
            setFIRA_ec(tti, {'sac_amp'}, sacs_(1,7));
%             rewsize = getFIRA_ec(tti, {'rew_off'}) - getFIRA_ec(tti, {'rew_on'});
%             setFIRA_ec(tti, {'reward'}, rewsize);
            setFIRA_ec(tti, {'sac_on_offline'}, getFIRA_ec(tti, sfixoff)+ sacs_(1,1));
            if nSacs>1 && ~isinf( sacs_(2,1) ) && ~isnan( sacs_(2,1) )
                t_eyeshift = sacs_(2,1) - sacs_(1,1);
                d_eyeshift = sqrt( (sacs_(2,5) - sacs_(1,5))*(sacs_(2,5) - sacs_(1,5)) + (sacs_(2,6) - sacs_(1,6))*(sacs_(2,6) - sacs_(1,6)) );
            end
            switch taskid
                case {2,3,4,5}  % memory & visual & ring
                    %Pick the active target to see if correct
                    ActiveTarg=getFIRA_ec(tti, {'trialid'});
                    if mod(ActiveTarg,2)==0
                        TargX=t1_x;
                        TargY=t1_y;
                        ActiveAngle=t1_angle;
                        InactiveAngle=t2_angle;

                    else
                        TargX=t2_x;
                        TargY=t2_y;
                        ActiveAngle=t2_angle;
                        InactiveAngle=t1_angle;


                    end
                    setFIRA_ec(tti, {'active_target'}, ActiveAngle);

                    delta = sqrt( (TargX - sacs_(1,5))*(TargX - sacs_(1,5)) + (TargY - sacs_(1,6))*(TargY - sacs_(1,6)) );
                    tarwin = sqrt(TargX*TargX + TargY*TargY)/1.2;
                    if delta<tarwin && isnan(getFIRA_ec(tti, {'OLerror'})) && ~isnan(sacs_(1,1)) && (abs(degAngDiff(ActiveAngle,mod(atan2(sacs_(1,6), sacs_(1,5))*180/pi+ 360,360)))<20) && (t_eyeshift>100 || d_eyeshift<=tarwin)
              
                        setFIRA_ec(tti, {'Offline_Score'}, 1);
                    elseif isnan(sacs_(1,1))
                        setFIRA_ec(tti, {'Offline_Score'}, 3);
                        
                    elseif ~isnan(getFIRA_ec(tti, {'OLerror'}))
                 
                        setFIRA_ec(tti, {'Offline_Score'}, 0);
                        
                        %Went to the sample, not hte target
                    elseif (abs(degAngDiff(InactiveAngle,sac_angle))<abs(degAngDiff(ActiveAngle,sac_angle))||(getFIRA_ec(tti, {'rew_on'}))==0)
                    setFIRA_ec(tti, {'Offline_Score'}, 2);
                    else 
                    setFIRA_ec(tti, {'Offline_Score'}, -1);

                    end
%               
            end

        end    
        currentscore = getFIRA_ec(tti,{'Offline_Score'});
        if ~isnan(getFIRA_ec(tti, {'OLcorrect'}))
            OLscore = 1;
% note: 
%   for recordings before 4/22/09, All_OffCD is the feedback for both correct and error trials. 
%     Reward follows All_OffCD in correct trials
%   for recordings afterward, 
%     in error trials, FDBKONCD is dropped 400 ms after TRGACQUIRECD, All_OffCD dropped 400 ms later.
%     in correct trials, REW DIO dropped 800 ms after TRGACQUIRECD
%             if (getFIRA_ec(tti, 'taskid')==3)
%                 setFIRA_ec(tti, {'fdbkon'}, getFIRA_ec(tti, {'All_Off'}));
%             end
        elseif ~isnan(getFIRA_ec(tti, {'OLerror'}))
            OLscore = 0;
        elseif ~isnan(getFIRA_ec(tti, {'OLncerr'}))
            OLscore = -1;
        elseif ~isnan(getFIRA_ec(tti, {'OLbrfix'}))
            OLscore = -2;
        else
            OLscore = -3;   %will this ever happen?
        end
        setFIRA_ec(tti, {'OLscore'}, OLscore);
        if currentscore == OLscore
            setFIRA_ec(tti, {'scorematch'}, 1);
        else
            setFIRA_ec(tti, {'scorematch'}, 0);
        end

    end

    
    
end
