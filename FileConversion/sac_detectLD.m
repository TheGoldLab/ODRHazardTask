function sac_ = sac_detectLD(h_eye, v_eye, store_rate, num_sacs)
% sac_detectLD.m   Decects the first saccade in a trace and returns other dynamics of the saccade too.
% function sac_ = sac_detectLD(h_eye, v_eye, store_rate, num_sacs)
%   v_threshold: velocity threshold (default: 30)  
%   a_threshold: acceleration threshold (default: 8000)
% This function detects a saccade using velocity and acceleration criterions
%   Since velocity and acceleration have been calucated, peak veloicty, acceleration and
%   deceleration are also returned.
%
% Note: This function requries the eye to be in double precision and MRDD returns the eye
%   as signle percision

%1/29/2018, Changed line 188 so that the selected saccade had a minumum
%duration and ended within 15 but farther than 5 units from the center,
%rather than having a minimum displacement difference from center, so that
%if the first sacade (or blink mascarading as a saccade)
%moved to the radius of the circle but the second
%corrected to the appropriate angle, that the second sacade would not be
%eliminated from consideration due to a small delta radius

% 8/21/2009, LD, added requirement of minmum sac duration to counter
% delta-like noise in eye traces.
% 3/08/2007, LD, renamed the return variables for compatibility with spm in
%       Josh's lab
% 12/16/03, LD. renamed file. use different threshold algorithms
% 12/9/03, LD. include vel and accel thresholds as input arguments
% modified by Long Ding 11/30/2003
% from Nicholas Port's program saccade_detect.m (created 9/24/99)
% 
%
% Changes, Nicholas Port, 3/17/98
% Was having problems with micro-saccades in the fixation period.  Now the displacement is calculated from
% 50 ms into the trial if there was a microsaccade

%Inits
sac_ = [];
sac_onset  = 0;
sac_offset = 0;
cross1     = [];
peak_vel   = [];
peak_accel = [];
peak_decel = [];
sac_sizeH   = [];
sac_sizeV   = [];
startH = [];
startV = [];
t_pv       = -1;
t_pa       = [];
t_pd       = [];
nSac = 0;
h_eye = h_eye';
v_eye = v_eye';

MIN_SACCADE_SIZE   = 2.5;
MAX_SACCADE_DURATION = 200;
MIN_SACCADE_DURATION = 10;
pre_dur = 100;
post_dur = 200;

% figure
% h = gcf;

[m n] = size( h_eye );

% check if the signal length is long enough
% filtfilt requires 3 times the filter length
if (n<150)       
    return;
end
h_baseline = mean(h_eye(1:10));
v_baseline = mean(v_eye(1:10));
h_eye = h_eye - h_baseline;
v_eye = v_eye - v_baseline;
% Calculate the displacement of the eye
t(1:n) = sqrt ( ( double(h_eye(1:n))).^2 + ...
                ( double(v_eye(1:n))).^2 );
h_eye = h_eye + h_baseline;
v_eye = v_eye + v_baseline;
% t(1:n) = abs(double(h_eye(1:n)));
            
% Setup the filter
a = 1;
b = fir1 ( 50, [0.30] );
% calculate velocity of the eye
vel(1:n)=gradient(t)*1000;  %deg/sec
% calculate acceleration of the eye
accel(1:n) = gradient(vel)*1000; %deg/sec/sec
% Filtering
vel = filtfilt (b, a, vel);
accel = filtfilt(b,a, accel);
% deteriming vel threshold
sd_vel = std(vel);
mean_vel = abs(mean(vel));
th_vel = 30;%mean_vel+2*sd_vel;
% th_vel = mean_vel+2*sd_vel;
% deteriming accel threshold
sd_accel = std(accel);
mean_accel = abs(mean(accel));
th_accel = 8000;
% th_accel = mean_accel+2*sd_accel;

% find candidate saccades based on vel and accel thresholds
% use the first crossing for sac
cand = find(vel>=th_vel & accel>=th_accel );
% remove crossings that are <50 ms from the end
j = find(cand<n-50);
cand = cand(j);
if isempty(cand)
    % No sac detected
    sac_ = [NaN NaN NaN NaN double(h_eye(end)) double(v_eye(end)) NaN NaN NaN NaN NaN];
    return
else
    dt_sec = double(1/store_rate);
    dt_ms = 1000*dt_sec;

    [m n_cand] = size(cand);
    % find potential sac onset. remove duplicates
    onset = zeros(1,n_cand);
    for i=1:n_cand
        vel_new = vel(1:cand(i));
        j = find(vel_new<=0);
        if ~isempty(j)
            % Sac onset is when vel cross zero before crossing, if searched
            % back in time
            m = length(j);
            onset(i) = j(m);
        end
    end
    [onset_unique, m_unique, n_unique] = unique(onset);  % remove duplicates
    cand = cand(m_unique);
    m = find(onset_unique>0);   % remove zero onset
    onset_unique = onset_unique(m);
    cand = cand(m);
    % Find the Peak Velocity, Peak Acceleration, and Peak Deceleration for each onset,
    % search between two onset times, or the last onset+ 200 ms(if not at the end of signal).
    m = length(onset_unique);
    if m<1
    % No sac detected
        return
    end
    onset_unique(m+1) = min( onset_unique(m)+MAX_SACCADE_DURATION, n );
    for i=1:m
        [peak_vel(i), t_pv(i)] = max(vel( onset_unique(i):onset_unique(i+1) ));
        [peak_accel(i), t_pa(i)] = max(accel( onset_unique(i):onset_unique(i+1) ));
%8-23-09        if onset_unique(i+1)>cand(i)+1
%8-23-09             [peak_decel(i), t_pd(i)] = min(vel( cand(i)+1:onset_unique(i+1) ));
        if onset_unique(i+1)>onset_unique(i)+1
            [peak_decel(i), t_pd(i)] = min(vel( onset_unique(i)+1:onset_unique(i+1) ));
        else
            peak_decel(i) = 0; t_pd(i)= 0; 
        end
        t_pv(i) = t_pv(i)+ onset_unique(i);
        t_pa(i) = t_pa(i)+ onset_unique(i);
        t_pd(i) = t_pd(i)+ onset_unique(i);
        % Find when the Saccade Ends, Starting from Peak Deceleration,
        % Look Forward in Velocity till it becomes below a constant
%8-23-09         j = find( vel( t_pd(i):onset_unique(i+1) ) < 30 );
        j = find( abs(vel( t_pv(i):onset_unique(i+1))) < 20 );
        if isempty(j)
            sac_offset(i) = -1;
            sac_size(i) = 0;
            sac_dur(i) = 0;
            sac_endDisp(i)=0;
        else
%8-23-09             i_goodj = find(j>20);
            i_goodj = find(j>10);
            if isempty(i_goodj)
                sac_offset(i) = -1;
                sac_size(i) = 0;
                sac_dur(i) = 0;
                sac_endDisp(i)=0;
                
            else
                j = j(i_goodj);
%                 sac_offset(i) = t_pd(i)+ j(1) - 1;
                sac_offset(i) = t_pv(i)+ j(1) - 1;
                sac_dur(i) = (sac_offset(i) - onset_unique(i))*dt_ms;
                sac_size(i) = t(sac_offset(i))-t(onset_unique(i));
                sac_endDisp(i)=abs(t(sac_offset(i)));
                startH(i) = double(h_eye(onset_unique(i)));
                startV(i) = double(v_eye(onset_unique(i)));
                sac_sizeH(i) = double(h_eye(sac_offset(i))) - double(h_eye(onset_unique(i)));
                sac_sizeV(i) = double(v_eye(sac_offset(i))) - double(v_eye(onset_unique(i)));
            end
        end
    end
    %CTRL Z to here to return to normal.  also Min dur was 10, now 25
    j = find(sac_endDisp<15 & sac_endDisp>5 & sac_dur>MIN_SACCADE_DURATION); %sac_size>MIN_SACCADE_SIZE 
    if isempty(j)
        sac_ = [NaN NaN NaN NaN double(h_eye(end)) double(v_eye(end)) NaN NaN NaN NaN NaN];
        return
    else
        nj = length(j);
        for i=1:nj
            rawsac(i) = sum(vel(onset_unique(j(i)):sac_offset(j(i)))*dt_sec);
        end
        dur = (sac_offset(j) - onset_unique(j))*dt_ms;
        sac_ = [onset_unique(j)*dt_ms; ...      % latency in ms
                dur; ...                     	% duration in ms
                peak_vel(j); ...                % peak velocity in deg/sec
                sac_size(j)./(dur*dt_sec); ...   % avg velocity in deg/sec
                sac_sizeH(j) + startH(j); ...   % end_x in deg
                sac_sizeV(j) + startV(j); ...   % end_y in deg
                rawsac; ...                     % point-by-point sum of sac trajectory (deg)
                sac_size(j); ...                % vector sac size (deg)
                sac_sizeH(j); ...               % sacsize_H (deg)
                sac_sizeV(j); ...               % sacsize_V (deg)
                t_pv(j)*dt_ms ...               % time of peak velocity (ms)
                ];                   
        [m,n]=size(sac_);
        if n>num_sacs
            sac_ = sac_(:, [1:num_sacs]);
        end           
        sac_ = sac_';
                
% 
%         sac_onset = onset_unique(j);
%         sac_offset = sac_offset(j);
%         sac_size = sac_size(j);
%         startH = startH(j);
%         startV = startV(j);
%         sac_sizeH = sac_sizeH(j);
%         sac_sizeV = sac_sizeV(j);
%         peak_vel = peak_vel(j);
%         peak_accel = peak_accel(j);
%         peak_decel = peak_decel(j);
%         t_pv = t_pv(j);
    end
end

% From findSaccadesA
% Returns an nx8 vectors, rows are saccades found, columns are:
%  1 .. latency
%  2 .. duration (ms)
%  3 .. maximum speed
%  4 .. average speed
%  5 .. end point x
%  6 .. end point y
%  7 .. distance raw  (the actual length of the saccade, calculated sample-by-sample)
%  8 .. distance vect (the length of the saccade calculated as the linear
%  vector)