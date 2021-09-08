function [ allData ] = KAS_Parser_794_DelayNeural( Cutoff, Type, data, theFile )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%Code to analyize the bias in responses based on SACP, current target, and
%sample mean

%Any .mat files you want to parse should be in the DataToParseFolder

%This one also registers whether or not the mean was shown during the
%decision period

%Gather Data
data_path = 'C:\Users\kyras\Desktop\UPenn 17-18\Matlab_AODR\ODRHazardTask-main\ConvertedFiles';%AVData/794'; % DatatoParse2';
% files = dir(fullfile(data_path, '*.mat'));
% nfiles = length(files);
cd(data_path);

% loop through the files and collect the following data matrix, 
%  rows are trials, columns are:
%  1. session
%  2. Trial index within session
%  3. Hazard rate
%  4. Noise (not currently used)
%  5. TACP(Actual)
%  6. T1_Angle
%  7. T2_Angle
%  8. Sample Angle
%  9. LLR for T1>T2
%  10. Offline Score: 0=err, 1=correct, 2=sample, -1=ncerr
%  11. Correct Target Shown? 0=no, 1=yes
%  12. Response_angle
%  13. Active_Angle 
%  14. Choice
%  15. Choice~=Correct(previous trial) aka switch
%  16. LLR for Change point
%  17. Online score
%  18. Correct 
%  19. Sample On time
%  20. Sample Off time
%  21. FP Off time
%  22. Sac On Time
%  23. Targ Acq Time
%  24. Number of units this session
%  25-end, spike times for units this session


%%
allData = [];
spikies=[];

% Please use the product of this script in conjunction with original .m data files
% with VisualizeTrial.m in ordr to examine individual trials

% Loop through the files
nfiles=1;
for currentFile = 1:nfiles

   % select "good" trials (no broken fixation) (Online score is not BF) and
   % is a AdODR task (taskid>=2);
   if Type==1
       task=data.ecodes.data(:,29)==1;
   else
       task=data.ecodes.data(:,29)>=2;
   end
   UnbrokenTrials = data.ecodes.data(:,41)>=-1 & task ;%& ~isnan(Choice);
   UnbrokenTrials_Index=find(UnbrokenTrials~=0);
   
   
      %Which target was picked?

   InactiveTarg=nan(size(data.ecodes.data(:,41)));
   InactiveTarg(data.ecodes.data(:,35)==data.ecodes.data(:,36))=data.ecodes.data(data.ecodes.data(:,35)==data.ecodes.data(:,36),37);
   InactiveTarg(data.ecodes.data(:,35)~=data.ecodes.data(:,36))=data.ecodes.data(data.ecodes.data(:,35)~=data.ecodes.data(:,36),36);


Choice2=nan(size(data.ecodes.data(UnbrokenTrials,41)));
  ActiveTarg=data.ecodes.data(UnbrokenTrials,35);
  RespAng=data.ecodes.data(UnbrokenTrials, [39]);
  Sample=data.ecodes.data(UnbrokenTrials, [38]);
  InactiveTarg=InactiveTarg(UnbrokenTrials);

for i=1:sum(UnbrokenTrials)
    if abs(degAngDiff(ActiveTarg(i),RespAng(i)))<=Cutoff
        Choice2(i)=ActiveTarg(i);
    elseif abs(degAngDiff(InactiveTarg(i),RespAng(i)))<=Cutoff
            Choice2(i)=InactiveTarg(i);

    end
end
   
test=find(UnbrokenTrials==1);
test2=~isnan(Choice2);
test3=test(test2);
temp=zeros(size(UnbrokenTrials));
temp(test3)=1;
UnbrokenTrials=(temp==1);
 UnbrokenTrials_Index=find(UnbrokenTrials~=0);
 
   Choice2=nans(size(data.ecodes.data(UnbrokenTrials,41)));
  ActiveTarg=data.ecodes.data(UnbrokenTrials,35);
  RespAng=data.ecodes.data(UnbrokenTrials, [39]);
  Sample=data.ecodes.data(UnbrokenTrials, [38]);
  
     InactiveTarg=nans(size(data.ecodes.data(:,41)));
   InactiveTarg(data.ecodes.data(:,35)==data.ecodes.data(:,36))=data.ecodes.data(data.ecodes.data(:,35)==data.ecodes.data(:,36),37);
   InactiveTarg(data.ecodes.data(:,35)~=data.ecodes.data(:,36))=data.ecodes.data(data.ecodes.data(:,35)~=data.ecodes.data(:,36),36);

  InactiveTarg=InactiveTarg(UnbrokenTrials);

for i=1:sum(UnbrokenTrials)
    if abs(degAngDiff(ActiveTarg(i),RespAng(i)))<=Cutoff
        Choice2(i)=ActiveTarg(i);
    elseif abs(degAngDiff(InactiveTarg(i),RespAng(i)))<=Cutoff
            Choice2(i)=InactiveTarg(i);

    end
end
   %Hazar Rate and Noise:
   Haz=data.ecodes.data(UnbrokenTrials, 32);
   Noise=data.ecodes.data(UnbrokenTrials, 33);
   
   
   %LLR for T1
   
      T1=max([normpdf(data.ecodes.data(UnbrokenTrials,38),data.ecodes.data(UnbrokenTrials,36),Noise)';normpdf(data.ecodes.data(UnbrokenTrials,38)+360,data.ecodes.data(UnbrokenTrials,36),Noise)';normpdf(data.ecodes.data(UnbrokenTrials,38)-360,data.ecodes.data(UnbrokenTrials,36),Noise)']); 
      T2=max([normpdf(data.ecodes.data(UnbrokenTrials,38),data.ecodes.data(UnbrokenTrials,37),Noise)';normpdf(data.ecodes.data(UnbrokenTrials,38)+360,data.ecodes.data(UnbrokenTrials,37),Noise)';normpdf(data.ecodes.data(UnbrokenTrials,38)-360,data.ecodes.data(UnbrokenTrials,37),Noise)']); 

   
   LLR=log(T1'./T2');
   
   %Were the choices visbible?
   Choice_Vis=~isnan(data.ecodes.data(UnbrokenTrials,8)) & (data.ecodes.data(UnbrokenTrials,8))<(data.ecodes.data(UnbrokenTrials,9));
   
  %Was the choice equal to the correct target of the previous trial?
  ActiveTarg=data.ecodes.data(UnbrokenTrials,35);
 
  Switch=(Choice2(2:end)~= ActiveTarg(1:end-1) & ~isnan(Choice2(2:end)));
  Switch=[nan;Switch];

  %LLR for change point
  InactiveTarg=nans(sum(UnbrokenTrials), 1);
  ActiveTarg=data.ecodes.data(UnbrokenTrials,35);
  Active1=data.ecodes.data(UnbrokenTrials,35)==data.ecodes.data(UnbrokenTrials,36);
  Targ1=data.ecodes.data(UnbrokenTrials,36);
  Targ2=data.ecodes.data(UnbrokenTrials,37);
  InactiveTarg(Active1)=Targ2(Active1);
    InactiveTarg(~Active1)=Targ1(~Active1);
    UnbrokenSamples= data.ecodes.data(UnbrokenTrials,38); 
    
   changeTarg=max([normpdf(UnbrokenSamples(2:end),InactiveTarg(1:end-1),Noise(1:end-1))';normpdf(UnbrokenSamples(2:end)+360,InactiveTarg(1:end-1),Noise(1:end-1))';normpdf(UnbrokenSamples(2:end)-360,InactiveTarg(1:end-1),Noise(1:end-1))']); 
   SameTarg=max([normpdf(UnbrokenSamples(2:end),ActiveTarg(1:end-1),Noise(1:end-1))';normpdf(UnbrokenSamples(2:end)+360,ActiveTarg(1:end-1),Noise(1:end-1))';normpdf(UnbrokenSamples(2:end)-360,ActiveTarg(1:end-1),Noise(1:end-1))']); 

    
  LLRCP=log(changeTarg./SameTarg)';
    LLRCP=[nan;LLRCP];
   
   
   % figure out trials since change points in unbroken trials
            %**** We are pretending broken trials straight up didn't happen
   % first extra-dimensional
   CPS = [0;find(diff(ActiveTarg)~=0);sum(UnbrokenTrials)];
   TACP = nans(sum(UnbrokenTrials),1);
   for ii = 1:length(CPS)-1
      inds = CPS(ii)+1:CPS(ii+1);
      TACP(inds) = 1:length(inds);
        %Allows you to assign a matrix to the slots of a matrix
   end


imptUnits=(mod(data.spikes.id,10)~=0);
%imptUnits=ones(size(data.spikes.id));
numNeurons=sum(imptUnits);
   if sum(UnbrokenTrials)>0
   % now collect it all   
   allData = cat(1, allData, [ ...
      currentFile.*ones(sum(UnbrokenTrials), 1), ...
      UnbrokenTrials_Index,...
      Haz,...
      Noise,...
      TACP, ...
      data.ecodes.data(UnbrokenTrials, [36 37 38]),... 
      LLR,...
      data.ecodes.data(UnbrokenTrials, [40]),...
      Choice_Vis,...
      data.ecodes.data(UnbrokenTrials, [39 35]), ... 
      Choice2,...
      Switch,...
      LLRCP,...
      data.ecodes.data(UnbrokenTrials,[41]),...
      Choice2==data.ecodes.data(UnbrokenTrials, [35]),...
      data.ecodes.data(UnbrokenTrials, [7 9 10 12,13]),...
      numNeurons.*ones(sum(UnbrokenTrials), 1)
      ]);   
  if ~isempty(data.spikes.data)
  spikies=[spikies;data.spikes.data(UnbrokenTrials,find(imptUnits==1))];
  else
      spikies=[spikies;cell(sum(UnbrokenTrials),1)];
  end
end
end


spikelist={};
neuronnames=data.spikes.id(find(imptUnits==1));
for i=1:numNeurons
    name=['spiketimes_',num2str(neuronnames(i))];
    spikelist=[spikelist,name];
end
if numNeurons==0
    spikelist={'none'};
end




Interm=num2cell(allData);
Interm=[Interm,spikies];
%
if ~isempty(allData)
for i=1:nfiles
Interm((allData(:,1)==i),1)={theFile};
end
% Make a Table

thenames={'Session','Trial_Within_Session','H_Rate',...
                                             'Noise', 'TACP_Actual', 'T1_Angle',...
                                             'T2_Angle', 'Sample_Angle', 'LLR', 'Offline_Score', 'Correct_On',...
                                             'Response_Angle', 'Active_Angle', 'Choice',...
                                             'Switch', 'LLR_for_Change', 'OnlineScore','Correct','SampleOn','SampleOff','FPOff','SacOn','TargAcq','Num_Neuron'};
thenames=[thenames,spikelist];

allData=cell2table(Interm,'VariableNames',thenames);
else
    allData=[];
end

end

