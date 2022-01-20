%Where your MAT Files are found
data_path = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\SortedConverted';
%'C:\Users\kyras\Desktop\UPenn 17-18\Matlab_AODR\ODRHazardTask-main\ConvertedFiles';%AVData/794'; % DatatoParse2';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

Collection=cell(nfiles,3);
close all

addpath 'C:\Users\alice\Documents\Projects\ODRHazardTask\Analysis\KyraPipeline20210129'
addpath 'C:\Users\alice\Documents\Projects\ODRHazardTask\FileConversion'



for currentFile = nfiles
    %(responses within this many degree of target are considered a response, rest are NaN)
    ResponseCutoff=30;
    theFile= files(currentFile).name;
    temp=load(theFile);
    data=temp.data;
    [~,b,~]=fileparts(theFile);
    mkdir(b)
    
    %Parse the files for Memory and Task trials, will gather all the
    %relevant information
    MemTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,1,data,theFile,data_path);
    TaskTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,2,data,theFile,data_path);% KAS_Parser_794_recoded(ResponseCutoff);
    if ~isempty(TaskTrials)
        numneuro=unique(TaskTrials.Num_Neuron);
    
    %Start making figures/doing analysis
    totalFigs=1+7*numneuro;
    cd(b)
    TACP_Cuttoff=5;
    
          [Fits,TrialNum,evInd]=FindHRateQuant(TaskTrials,1+totalFigs*(currentFile-1),b);

    % Subsets are 1= Corr/Err, 2=Stay/Switch, 3=Chose 1/Chose 2
    %AFR and STDFR contain the mean/StD FR across all the trials for a
    %given selection, whereas PeriodMeanFR/PeriodSTDFR give the mean/std
    %FR across a given period (1=vis, 2=mem, 3=motor, 4=rew).
    %NumTriBreakdown contains the number of trials that went into the
    %averages
    %AvgFR organized by (Ev, Hrate,Active,SubsetSelection(eg corr/err),neuron number)
    Alice=1; %If alice is using this/wants the bar graphs from the peridod summary
    for subset=1:3
       [AvgFr{subset},STDFR{subset}, NumTriBreakdown{subset}]=MonkeySpikeAnalysisSubset(TaskTrials,subset,(numneuro*(subset))+2+totalFigs*(currentFile-1),b,TACP_Cuttoff);  
        [PeriodMeanFR{subset},PeriodStdFR{subset}]=MonkeySpikeAnalysisSubset_Period_Summary(TaskTrials,subset,(numneuro*(subset))+5+totalFigs*(currentFile-1),b,TACP_Cuttoff,Alice);
        
    end
    if ~isempty(MemTrials)
         [directionFits_MGS, MFR_MGS, SFR_MGS, NumTri_MGS]=MonkeySpikeDirectionality(MemTrials,numneuro*(subset+1)+8+totalFigs*(currentFile-1),b);
    end
    
    
    [AvgFRTotal_LLRCh, STDFRTotal_LLRCh, NumTri_LLRCh ] =MonkeySpikeAnalysisSubset2(TaskTrials);
    end

    
    
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
    close all
    
end

