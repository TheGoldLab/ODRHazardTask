% function fDebugFIRA
%  For lokking at FIRA while debugging the programs
%  
%  01/23/12 Taka
%--------------------------------------------------------------------------
function fDebugFIRA
global FIRA

correct       = FIRA.ecodes.data(:,find(strcmp(FIRA.ecodes.name,'correct')));
OLscore       = FIRA.ecodes.data(:,find(strcmp(FIRA.ecodes.name,'OLscore')));
dot_coh       = FIRA.ecodes.data(:,find(strcmp(FIRA.ecodes.name,'dot_coh')));
tmp_Idx.OLcor = 17; tmp_Idx.OLerr = 18; tmp_Idx.OLnce = 19; tmp_Idx.OLbrf = 20;
OLtime        = FIRA.ecodes.data(:,[tmp_Idx.OLcor tmp_Idx.OLerr tmp_Idx.OLnce tmp_Idx.OLbrf]);
estim         = FIRA.ecodes.data(:,find(strcmp(FIRA.ecodes.name,'estim_type')));

correct_ = correct(estim==1);
OLscore_ = OLscore(estim==1);
OLtime_  = OLtime(estim==1,:);
dot_coh_ = dot_coh(estim==1);

if sum(OLscore_ == correct_)<length(OLscore_)
    disp([mfilename, ': the mismatch between OLscore and correct occures.']);
%     keyboard;
end

keyboard;