data_path = '/Users/KAS/Documents/MATLAB/ODR Task Analysis/AdaptODRMonkey/DatatoParse2';%AVData/794'; % DatatoParse2';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

Collection=cell(nfiles,3);
close all
Lapses=[],
for currentFile = 34% 1:nfiles
ResponseCutoff=30; %(responses within this many degree of target are considered a response, rest are NaN)
theFile= files(currentFile).name;
temp=load(theFile);
data=temp.data;
[~,b,~]=fileparts(theFile);
mkdir(b)



MemTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,1,data,theFile);

TaskTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,2,data,theFile);% KAS_Parser_794_recoded(ResponseCutoff);
  numneuro=unique(TaskTrials.Num_Neuron);

totalFigs=1+5*numneuro;
cd(b)
  [Fits,TrialNum]=FindHRateQuant(TaskTrials,1+totalFigs*(currentFile-1));
    saveas(gcf,[b,'_PsychCurve.pdf'])
    if istable(Fits)
    Lapses(currentFile)=Fits.Fit_Lapse(1);
   save([b,'_Fits'],'Fits');
    end

%  MonkeySpikeAnalysisTotal(TaskTrials,2+totalFigs*(currentFile-1),b);
%   
% for subset=1:3
%  NumTriBreakdown{subset}=MonkeySpikeAnalysisSubset(TaskTrials,subset,(numneuro*(subset))+2+totalFigs*(currentFile-1),b);
% end
% if ~isempty(MemTrials)
% [directionFits,NumTri]=MonkeySpikeDirectionality(MemTrials,numneuro*(subset+1)+2+totalFigs*(currentFile-1),b)
% end
% 
% Collection{1,1}=theFile;
% Collection{1,2}=MemTrials;
% Collection{1,3}=TaskTrials;
% Collection{1,4}=Fits;
% Collection{1,5}=TrialNum;
% Collection{1,6}=NumTriBreakdown{1};
% Collection{1,7}=NumTriBreakdown{2};
% Collection{1,8}=NumTriBreakdown{3};
% if ~isempty(MemTrials)
% Collection{1,9}=directionFits;
% Collection{1,10}=NumTri;
% 
% 
% Collection=cell2table(Collection,'VariableNames',{'Name','Mem_Trials','Task_Trials','Psycho_fit', 'Trial_Nums','TriNum_Corr','TriNum_Switch','TriNumChoice','PrefDir_Fits','Dir_Trials'});
% save([b,'_Parsed'],'Collection');
 cd(data_path);
% clear Collection
% 
% else
%    Collection=cell2table(Collection,'VariableNames',{'Name','Mem_Trials','Task_Trials','Psycho_fit', 'Trial_Nums','TriNum_Corr','TriNum_Switch','TriNumChoice'});
% save([b,'_Parsed'],'Collection');
% cd(data_path);
% clear Collection 
% end
end

