%Where your MAT Files are found
data_path = 'C:\Users\kyras\Desktop\UPenn 17-18\Matlab_AODR\ODRHazardTask-main\ConvertedFiles';%AVData/794'; % DatatoParse2';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);
AvgFr={};
STDFR={};
NumTriBreakdown={};
Collection=cell(nfiles,3);
close all
pvalsGS=[];
CellNum=1;
%MegaDataBase=table;
%MegaDataBase.Properties.VariableNames={'Monkey','DirV','DirV_p','DirMem','DirMem_p','DirMot','DirMot_p','DirR','DirR_p','Ev','Ev_p','HazAct','HazAct_p','Delta_T1','Delta_N','Delta_T2','FitH_L','FitH_N','Bias_L','Bias_N'};
%%
TaskTrialsOfInt=[];
for currentFile =  4:nfiles
    %nfiles
    Collection=cell(1,3);
    %(responses within this many degree of target are considered a response, rest are NaN)
    ResponseCutoff=30;
    theFile= files(currentFile).name;
    temp=load(theFile);
    data=temp.data;
    [~,b,~]=fileparts(theFile);
    cd(data_path);

    mkdir(b)
    
    %Parse the files for Memory and Task trials, will gather all the
    %relevant information
    MemTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,1,data,theFile,data_path);
    TaskTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,2,data,theFile,data_path);% KAS_Parser_794_recoded(ResponseCutoff);
   
    TaskTrialsOfInt=[TaskTrialsOfInt;TaskTrials(:,1:23)];
    
    if ~isempty(TaskTrials)
    numneuro=unique(TaskTrials.Num_Neuron);
    else
        numneuro=unique(MemTrials.Num_Neuron);
    end
    %Start making figures/doing analysis
    totalFigs=1+7*numneuro;
    cd([data_path,'\',b])
    TACP_Cuttoff=5;
    if ~isempty(TaskTrials)
            [Fits,TrialNum]=FindHRateQuant(TaskTrials,1+totalFigs,b,1,[data_path,'\',b]);

    % Subsets are 1= Corr/Err, 2=Stay/Switch, 3=Chose 1/Chose 2
    %AFR and STDFR contain the mean/StD FR across all the trials for a
    %given selection, whereas PeriodMeanFR/PeriodSTDFR give the mean/std
    %FR across a given period (1=vis, 2=mem, 3=motor, 4=rew).
    %NumTriBreakdown contains the number of trials that went into the
    %averages
    %AvgFR organized by (Ev, Hrate,Active,SubsetSelection(eg corr/err),neuron number)
    Alice=0; %If alice is using this/wants the bar graphs from the peridod summary
    if numneuro>0
    for subset=1%:3
       [a,q,c]=MonkeySpikeAnalysisSubset(TaskTrials,subset,(numneuro*(subset))+2+totalFigs,b,TACP_Cuttoff);  
        [PeriodMeanFR{subset,currentFile},PeriodStdFR{subset},HrateEffect{currentFile},ComparisonNeutral,SignRankLowInt,SignRankLowSlope,SignRankNInt,SignRankNSlope,IntRealLowDiff,IntRealNDiff,SlpeRealLowDiff,SlpeRealNDiff,pval_H_Center_Diff]=MonkeySpikeAnalysisSubset_Period_Summary(TaskTrials,subset,(numneuro*(subset))+5+totalFigs,b,TACP_Cuttoff,Alice);
        AvgFR{subset}=a;
        STDFR{subset}=q;
        NumTriBreakdown{subset}=c;
        %% THIS IS WHERE YOU LOOK AT H RATE OVEr ALL, STILL NEED TO HAVE DIRECTIONALITY CRITERIA
        %
        for Neuro=1:numneuro
for P=1:5

                Diff{currentFile,Neuro}(P)=abs(PeriodMeanFR{1,currentFile}{1,Neuro}(5,1,1,1,P)-PeriodMeanFR{1,currentFile}{1,Neuro}(5,1,2,1,P))-abs(PeriodMeanFR{1,currentFile}{1,Neuro}(5,2,1,1,P)-PeriodMeanFR{1,currentFile}{1,Neuro}(5,2,2,1,P));
            LowDiff{currentFile,Neuro}(P)=abs(PeriodMeanFR{1,currentFile}{1,Neuro}(5,1,1,1,P)-PeriodMeanFR{1,currentFile}{1,Neuro}(5,1,2,1,P));
            NDiff{currentFile,Neuro}(P)=abs(PeriodMeanFR{1,currentFile}{1,Neuro}(5,2,1,1,P)-PeriodMeanFR{1,currentFile}{1,Neuro}(5,2,2,1,P));
end
        end
    end
    
    if ~isempty(MemTrials) 
         [directionFits_MGS, MFR_MGS, SFR_MGS, NumTri_MGS, pval,Index,Comparison]=MonkeySpikeDirectionality(MemTrials,numneuro*(subset+1)+8+totalFigs,b);
    else
        pvals=repmat(NaN, numneuro, 4);
    end
   
    
for cells=1:numneuro
    %General Info
    MegaDataBase_Neural{CellNum,1}={theFile(1:2)}; %Monkey
    MegaDataBase_Neural{CellNum,2}={theFile(4:13)}; %Date
    MegaDataBase_Neural{CellNum,3}=cells;          %Cell in the file
    MegaDataBase_Neural{CellNum,4}={theFile(15:18)}; %Depth (not for all)
    %Directionality for each time period (Vis,Mem,Mot,Rew), the pval, and the index
    MegaDataBase_Neural{CellNum,5}=directionFits_MGS{1,cells}{1,1}{1,2};       %Fit Direction
    MegaDataBase_Neural{CellNum,6}=pval(cells,1);                              %pval of whether one dir has sig increased FR
    MegaDataBase_Neural{CellNum,7}=Index(cells,1);                             %index of Fit strength 
    MegaDataBase_Neural{CellNum,8}=directionFits_MGS{1,cells}{1,1}{2,2};
    MegaDataBase_Neural{CellNum,9}=pval(cells,2);
    MegaDataBase_Neural{CellNum,10}=Index(cells,2);
    MegaDataBase_Neural{CellNum,11}=directionFits_MGS{1,cells}{1,1}{3,2};
    MegaDataBase_Neural{CellNum,12}=pval(cells,3);
    MegaDataBase_Neural{CellNum,13}=Index(cells,3);
    MegaDataBase_Neural{CellNum,14}=directionFits_MGS{1,cells}{1,1}{4,2};
    MegaDataBase_Neural{CellNum,15}=pval(cells,4);
    MegaDataBase_Neural{CellNum,16}=Index(cells,4);
    %Effect of Evidence in aODR
    MegaDataBase_Neural{CellNum,17}=HrateEffect{1,currentFile}{cells,1}{2,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,18}=HrateEffect{1,currentFile}{cells,1}{2,4}; %Significance
    MegaDataBase_Neural{CellNum,19}=HrateEffect{1,currentFile}{cells,2}{2,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,20}=HrateEffect{1,currentFile}{cells,2}{2,4}; %Significance
    MegaDataBase_Neural{CellNum,21}=HrateEffect{1,currentFile}{cells,3}{2,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,22}=HrateEffect{1,currentFile}{cells,3}{2,4}; %Significance
    MegaDataBase_Neural{CellNum,23}=HrateEffect{1,currentFile}{cells,4}{2,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,24}=HrateEffect{1,currentFile}{cells,4}{2,4}; %Significance
    MegaDataBase_Neural{CellNum,25}=HrateEffect{1,currentFile}{cells,5}{2,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,26}=HrateEffect{1,currentFile}{cells,5}{2,4}; %Significance
   %Effect of H*Active in aODR
    MegaDataBase_Neural{CellNum,27}=HrateEffect{1,currentFile}{cells,1}{5,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,28}=HrateEffect{1,currentFile}{cells,1}{5,4}; %Significance
    MegaDataBase_Neural{CellNum,29}=HrateEffect{1,currentFile}{cells,2}{5,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,30}=HrateEffect{1,currentFile}{cells,2}{5,4}; %Significance
    MegaDataBase_Neural{CellNum,31}=HrateEffect{1,currentFile}{cells,3}{5,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,32}=HrateEffect{1,currentFile}{cells,3}{5,4}; %Significance
    MegaDataBase_Neural{CellNum,33}=HrateEffect{1,currentFile}{cells,4}{5,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,34}=HrateEffect{1,currentFile}{cells,4}{5,4}; %Significance
    MegaDataBase_Neural{CellNum,35}=HrateEffect{1,currentFile}{cells,5}{5,1}; %Effect of Evidence location
    MegaDataBase_Neural{CellNum,36}=HrateEffect{1,currentFile}{cells,5}{5,4}; %Significance
    % Fit H'rates and biases
   MegaDataBase_Neural{CellNum,37}=Fits.Fit_H(1); %Hrate_Low
   MegaDataBase_Neural{CellNum,38}=Fits.Fit_H(2); %Hrate_N
   MegaDataBase_Neural{CellNum,39}=Fits.Fit_H(2)-Fits.Fit_H(1); %Hrate_Shift
   MegaDataBase_Neural{CellNum,40}=Fits.LogBias(1); %BiasLow
   MegaDataBase_Neural{CellNum,41}=Fits.LogBias(2); %Bias_N
   MegaDataBase_Neural{CellNum,42}=Fits.LogBias(1)-Fits.LogBias(2); %Shift in bias (want to be POS)
   
   %SigDifferences In activity of T3 vs V/MGS
   %ComparionsNeutral{cell,Time,Active} %First row numbrs, second row location
   %Comparison{cell,Time}%First row numbers, second row location
   for time=1:3
       for loc=1:3
           [~,indexing]=min(abs(ComparisonNeutral{cells,time,1}{2,loc}-[Comparison{cells,time}{2,:}]));
           NeutralHolder=[ComparisonNeutral{cells,time,1}{1,loc};ComparisonNeutral{cells,2,2}{1,loc}];
           VMGSHolder=Comparison{cells,time}{1,indexing};
           [~,pvalLocDiff(time,loc)]=ttest2(NeutralHolder,VMGSHolder);
           NeutralMFR(time,loc)=nanmean(NeutralHolder);
           GS_MFR(time,loc)=nanmean(VMGSHolder);
       end
   end
      MegaDataBase_Neural{CellNum,43}=pvalLocDiff(1,1); %In Vis Period is activity different for T1?
      MegaDataBase_Neural{CellNum,44}=pvalLocDiff(1,2); %In Vis Period is activity different for Neutral?
      MegaDataBase_Neural{CellNum,45}=pvalLocDiff(1,3); %In Vis Period is activity different for T2?
      MegaDataBase_Neural{CellNum,46}=pvalLocDiff(2,1); %In Mem Period is activity different for T1?
      MegaDataBase_Neural{CellNum,47}=pvalLocDiff(2,2); %In Mem Period is activity different for Neutral?
      MegaDataBase_Neural{CellNum,48}=pvalLocDiff(2,3); %In Mem Period is activity different for T2?
      MegaDataBase_Neural{CellNum,49}=pvalLocDiff(3,1); %In Mot Period is activity different for T1?
      MegaDataBase_Neural{CellNum,50}=pvalLocDiff(3,2); %In Mot Period is activity different for Neutral?
      MegaDataBase_Neural{CellNum,51}=pvalLocDiff(3,3); %In Mot Period is activity different for T2?
      %Differences in Low/neutralH by Period
      MegaDataBase_Neural{CellNum,52}=Diff{currentFile,cells}(1); %In V period, what is diff of Low to neutral H (Positive indicates Low H has larger diff)
MegaDataBase_Neural{CellNum,53}=Diff{currentFile,cells}(2);
MegaDataBase_Neural{CellNum,54}=Diff{currentFile,cells}(3);
MegaDataBase_Neural{CellNum,55}=Diff{currentFile,cells}(4);
MegaDataBase_Neural{CellNum,56}=Diff{currentFile,cells}(5);
MegaDataBase_Neural{CellNum,57}=SignRankLowInt(cells,2);
MegaDataBase_Neural{CellNum,58}=SignRankLowSlope(cells,2);
MegaDataBase_Neural{CellNum,59}=SignRankNInt(cells,2);
MegaDataBase_Neural{CellNum,60}=SignRankNSlope(cells,2);
MegaDataBase_Neural{CellNum,61}=IntRealLowDiff(cells,2);
MegaDataBase_Neural{CellNum,62}=SlpeRealLowDiff(cells,2);
MegaDataBase_Neural{CellNum,63}=IntRealNDiff(cells,2);
MegaDataBase_Neural{CellNum,64}=SlpeRealNDiff(cells,2);
 MegaDataBase_Neural{CellNum,65}=Fits.LogLaps(1);
       MegaDataBase_Neural{CellNum,66}=GS_MFR(2,1); %In Mem Period is activity different for T1?
      MegaDataBase_Neural{CellNum,67}=GS_MFR(2,2); %In Mem Period is activity different for Neutral?
      MegaDataBase_Neural{CellNum,68}=GS_MFR(2,3);
       MegaDataBase_Neural{CellNum,69}=NeutralMFR(2,1); %In Mem Period is activity different for T1?
      MegaDataBase_Neural{CellNum,70}=NeutralMFR(2,2); %In Mem Period is activity different for Neutral?
      MegaDataBase_Neural{CellNum,71}=NeutralMFR(2,3);
            MegaDataBase_Neural{CellNum,72}=LowDiff{currentFile,cells}(2);
      MegaDataBase_Neural{CellNum,73}=NDiff{currentFile,cells}(2);
MegaDataBase_Neural{CellNum,74}=pval_H_Center_Diff(cells,1);
MegaDataBase_Neural{CellNum,75}=pval_H_Center_Diff(cells,2); 
    
    CellNum=CellNum+1;
end
    %
     Collection{1,1}=theFile;
     Collection{1,2}=MemTrials;
     Collection{1,3}=TaskTrials;
%      Collection{1,4}=Fits;
   %  Collection{1,5}=TrialNum;
    % Collection{1,6}=NumTriBreakdown{1};
    % Collection{1,7}=NumTriBreakdown{2};
    % Collection{1,8}=NumTriBreakdown{3};
    % if ~isempty(MemTrials)
    % Collection{1,9}=directionFits;
    % Collection{1,10}=NumTri;
    %
    %
%      Collection=cell2table(Collection,'VariableNames',{'Name','Mem_Trials','Task_Trials','Psycho_fit'});%, 'Trial_Nums','TriNum_Corr','TriNum_Switch','TriNumChoice','PrefDir_Fits','Dir_Trials'});
%      save([b,'_Parsed'],'Collection');
    cd(data_path);
%     % clear Collection
%     lowH_Haz(currentFile)=Fits{1,1};
%     lowH_Log(currentFile)=Fits{1,5};
%     highH_Haz(currentFile)=Fits{2,1};
%     highH_Log(currentFile)=Fits{2,5};
 %   pvalsGS=[pvalsGS;pval];
 FitDirs{1,currentFile}=pval;
    FitDirs{2,currentFile}=directionFits_MGS;
    FitDirs{3,currentFile}=Index;
    close all
    end
    end
end
%%
MegaDataBase_Neural.Properties.VariableNames={'Monkey','Date','Cell','Depth','DirFit_V','DirDiff_p_V','Index_V','DirFit_Mem','DirDiff_p_Mem','Index_Mem','DirFit_Mot','DirDiff_p_Mot','Index_Mot','DirFit_Rew','DirDiff_p_Rew','Index_Rew',...
    'Ev_V','Ev_p_V', 'Ev_Mem', 'Ev_p_Mem', 'Ev_Mot', 'Ev_p_Mot', 'Ev_Rew', 'Ev_p_Rew', 'Ev_Post', 'Ev_p_Post',...
    'HAct_V','HAct_p_V', 'HAct_Mem', 'HAct_p_Mem', 'HAct_Mot', 'HAct_p_Mot', 'HAct_Rew', 'HAct_p_Rew',  'HAct_Post', 'HAct_p_Post',...
    'HRate_L','HRate_N','HRate_Diff', 'Bias_L','Bias_N','Bias_Shift_POS',...
    'T1_Diff_V','N_Diff_V','T2_Diff_V','T1_Diff_Mem','N_Diff_Mem','T2_Diff_Mem','T1_Diff_Mot','N_Diff_Mot','T2_Diff_Mot'...
    'DiffH_V','DiffH_Mem','DiffH_Mot','DiffH_Rew','DiffH_Post'...
    'DiffInt_Low_Mem_p','DiffSlope_Low_Mem_p','DiffInt_N_Mem_p','DiffSlope_N_Mem_p'...
    'DiffInt_Low_Mem','DiffSlope_Low_Mem','DiffInt_N_Mem','DiffSlope_N_Mem', 'Log_Lapse'...
    'T1_GSFR_M','N_GSFR_M','T2_GSFR_M','T1_NeutralFR_M','N_NeutralFR_M','T2_NeutralFR_M','LowDiff','NDiff','Pval_Center_Low','Pval_Center_Neutral'};
% %% Summary Plots Population
% figure(1)
% scatter(MegaDataBase.DirDiff_p_Mem,MegaDataBase.Index_Mem)
% xlabel({'P-val of Directionality test, Mem'})
% ylabel('Directionality Index, Memory')
% title('Comparison of Directionality index to Directionaity Significance test')

%% Behavior
%%
figure(1)
subplot(1,2,2)
[a,b]=unique(MegaDataBase_Neural.HRate_Diff(1:77));
sortedb=sort(b(~isnan(a)));
plot(MegaDataBase_Neural.HRate_Diff(sortedb),'-ok','Linewidth',2)
hold on
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[0,0],'--k')
ylabel({'{\Delta}Subjective H rate'})
xlabel('Session Number')
title('H rate Difference over time')
subplot(1,2,1)
plot(MegaDataBase_Neural.Bias_Shift_POS(sortedb),'-ok','Linewidth',2)
hold on
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[0,0],'--k')
ylabel({'{\Delta}Bias'})
title('Bias Difference over time')
xlabel('Session Number')
set(gcf, 'Units', 'Inches')
 %set(gcf,'Position',[7.8000    1.9167    6.8583    6.2250]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
%saveas(gcf,[b,'_PsychometricFxn.pdf'])  
hold=MegaDataBase_Neural.Bias_Shift_POS(sortedb);
signrank(hold)
signrank(MegaDataBase_Neural.HRate_Diff(sortedb))
%%
figure(2)
plot(MegaDataBase_Neural.Log_Lapse(sortedb),'-ok','Linewidth',2)
hold on
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[0,0],'--k')
xlabel({'Session Number'})
ylabel('Lapse Rate')
nanmean(MegaDataBase_Neural.Log_Lapse(sortedb))
%%
figure (3)

subplot(1,2,1)
sortedb=sort(b);
plot(MegaDataBase_Neural.HRate_L(sortedb),'-or','Linewidth',2)
hold on
plot(MegaDataBase_Neural.HRate_N(sortedb),'-ob','Linewidth',2)
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[.5,.5],'--b')
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[.05,.05],'--r')
ylabel({'Subjective H rate'})
xlabel('Session Number')
title('H rate over time')
subplot(1,2,2)
plot(MegaDataBase_Neural.Bias_L(sortedb),'-or','Linewidth',2)
hold on
plot(MegaDataBase_Neural.Bias_N(sortedb),'-ob','Linewidth',2)
hold on
plot([0,length(MegaDataBase_Neural.HRate_Diff(sortedb))],[0,0],'--k')
ylabel({'Psychometric Curve Bias';})
title('Bias over time')
xlabel('Session Number')
legend({'Low H', 'Neutral H'}, 'Location', 'northeast')


%% Figure of LL distributions
Pcs=[0,5,10,10,15,15,20,15,10];
Pc2=[10,15,20,15,15,10,10,5,0];
Pms=[0,0,0,10,20,30,15,15,10];
Pm2=[10,15,15,30,20,10,0,0,0];
figure(5)
subplot(1,2,1)
 plot([1:9],Pcs,'-or', 'Linewidth',2)
hold on
plot([1:9],Pc2,'-ob','Linewidth',2)
ylabel('Prob of Cue given correct')
xlabel('Cue location')
legend({'9 correct', '1 correct'},'Location','south')
axis([1,9,0,30])
title('Ci Distributions')
subplot(1,2,2)
 plot([1:9],Pms,'-or', 'Linewidth',2)
hold on
plot([1:9],Pm2,'-ob','Linewidth',2)
ylabel('Prob of Cue given correct')
xlabel('Cue location')
legend({'9 correct', '1 correct'},'Location','south')
axis([1,9,0,30])
title('MM Distributions')


set(gcf, 'Units', 'Inches')
 set(gcf,'Position',[5.0729    3.5521   10.0854    4.3750]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,['Distributions.pdf'])  
%% Index colored by pvals
figure(6)
subplot(1,3,1)
hold on
if sum((MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05))>0
plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Neither Sig')
end
if sum(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05)>0
 plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5,'DisplayName','Mem Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05))>0
          plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),'o','MarkerEdgeColor','r',...
              'MarkerFaceColor','r',...
              'LineWidth',1.5, 'DisplayName', 'Vis Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05))>0
          plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),'o','MarkerEdgeColor',[0.6941    .5    0.8510],...
              'MarkerFaceColor',[0.6941    .5    0.8510],...
              'LineWidth',1.5, 'DisplayName', 'Vis and Mem Sig')
end
legend('AutoUpdate','off','Location','best')
xlabel({'Directionality Index, Vis'})
ylabel('Directionality Index, Memory')
plot([0,1],[0,1],'--k')

subplot(1,3,2)
hold on
if sum((MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05))>0
plot(MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Neither Sig')
end
if sum(MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05)>0
 plot(MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_Mot>.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5,'DisplayName','Memory Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05))>0
          plot(MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem>.05),'o','MarkerEdgeColor','c',...
              'MarkerFaceColor','c',...
              'LineWidth',1.5, 'DisplayName', 'Motor Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05))>0
          plot(MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),MegaDataBase_Neural.Index_Mem(MegaDataBase_Neural.DirDiff_p_Mot<.05 & MegaDataBase_Neural.DirDiff_p_Mem<.05),'o','MarkerEdgeColor',[0.6941    .5    0.8510],...
              'MarkerFaceColor',[0.6941    .5    0.8510],...
              'LineWidth',1.5, 'DisplayName', 'Mem and Motor Sig')
end
legend('AutoUpdate','off','Location','best')
xlabel({'Directionality Index, Mot'})
ylabel('Directionality Index, Memory')
plot([0,1],[0,1],'--k')


subplot(1,3,3)
hold on
if sum((MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05))>0
plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05),MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Neither Sig')
end
if sum(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05)>0
 plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05),MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_V>.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05),'o','MarkerEdgeColor','c',...
              'MarkerFaceColor','c',...
              'LineWidth',1.5,'DisplayName','Motor Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05))>0
          plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05),MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot>.05),'o','MarkerEdgeColor','r',...
              'MarkerFaceColor','r',...
              'LineWidth',1.5, 'DisplayName', 'Vis Sig')
end
if sum((MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05))>0
          plot(MegaDataBase_Neural.Index_V(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05),MegaDataBase_Neural.Index_Mot(MegaDataBase_Neural.DirDiff_p_V<.05 & MegaDataBase_Neural.DirDiff_p_Mot<.05),'o','MarkerEdgeColor',[0.6941    .5    0.8510],...
              'MarkerFaceColor',[0.6941    .5    0.8510],...
              'LineWidth',1.5, 'DisplayName', 'Vis and Motor Sig')
end
legend('AutoUpdate','off','Location','best')
xlabel({'Directionality Index, Vis'})
ylabel('Directionality Index, Motor')
plot([0,1],[0,1],'--k')

set(gcf, 'Units', 'Inches')
set(gcf,'Position',[0.0167    3.5521   15.4833    4.3750]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,['TuningIndex.pdf'])  


%% Low differences in Int/Slope vs Neutal differences in Int/Slope
figure(7)
subplot(1,2,1)
hold on
if sum((MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05))>0
plot(abs(MegaDataBase_Neural.DiffInt_Low_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05)),abs(MegaDataBase_Neural.DiffInt_N_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05)),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Neither Sig')
end
if sum((MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05))>0
plot(abs(MegaDataBase_Neural.DiffInt_Low_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05)),abs(MegaDataBase_Neural.DiffInt_N_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p>.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05)),'o','MarkerEdgeColor','r',...
              'MarkerFaceColor','r',...
              'LineWidth',1.5, 'DisplayName', 'Neutral Sig')
end
if sum((MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05))>0
disp(sum((MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05)))
    plot(abs(MegaDataBase_Neural.DiffInt_Low_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05)),abs(MegaDataBase_Neural.DiffInt_N_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05)),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5, 'DisplayName', 'Low Sig**')
          
        find((MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p>.05))
end
if sum((MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05))>0
disp(sum((MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05)))
    plot(abs(MegaDataBase_Neural.DiffInt_Low_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05)),abs(MegaDataBase_Neural.DiffInt_N_Mem(MegaDataBase_Neural.DiffInt_Low_Mem_p<.05 & MegaDataBase_Neural.DiffInt_N_Mem_p<.05)),'o','MarkerEdgeColor',[0.6941    .5    0.8510],...
              'MarkerFaceColor',[0.6941    .5    0.8510],...
              'LineWidth',1.5, 'DisplayName', 'Both Sig')
end
legend('Autoupdate','off','Location', 'best')
axis([0,30,0,30])
title('Diff by Active of Intercept for FR=m*Ev+b')
xlabel( 'Differnce in Active Intercepts: Low');
ylabel('Difference in Active Intercepts: Neutral');
plot([0,30],[0,30],'--k')

subplot(1,2,2)
hold on
if sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05))>0
plot(abs(MegaDataBase_Neural.DiffSlope_Low_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05)),abs(MegaDataBase_Neural.DiffSlope_N_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05)),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Low: p>.05; Neutral: p>.05')
end
if sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05))>0
plot(abs(MegaDataBase_Neural.DiffSlope_Low_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05)),abs(MegaDataBase_Neural.DiffSlope_N_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p>.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05)),'o','MarkerEdgeColor','r',...
              'MarkerFaceColor','r',...
              'LineWidth',1.5, 'DisplayName', 'Low: p>.05; Neutral: p<.05')
end
if sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05))>0
disp(sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05)))
    plot(abs(MegaDataBase_Neural.DiffSlope_Low_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05)),abs(MegaDataBase_Neural.DiffSlope_N_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p>.05)),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5, 'DisplayName', 'Low: p<.05; Neutral: p>.05**')
end
if sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05))>0
disp(sum((MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05)))
    plot(abs(MegaDataBase_Neural.DiffSlope_Low_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05)),abs(MegaDataBase_Neural.DiffSlope_N_Mem(MegaDataBase_Neural.DiffSlope_Low_Mem_p<.05 & MegaDataBase_Neural.DiffSlope_N_Mem_p<.05)),'o','MarkerEdgeColor',[0.6941    .5    0.8510],...
              'MarkerFaceColor',[0.6941    .5    0.8510],...
              'LineWidth',1.5, 'DisplayName', 'Low: p<.05; Neutral: p<.05')
end

Slope_SigDiff=signrank(abs(MegaDataBase_Neural.DiffSlope_Low_Mem)-abs(MegaDataBase_Neural.DiffSlope_N_Mem))
Int_SigDiff=signrank(abs(MegaDataBase_Neural.DiffInt_Low_Mem)-abs(MegaDataBase_Neural.DiffInt_N_Mem))

plot([0,.2],[0,.2],'--k')
% legend('Autoupdate','off','Location', 'best')
axis([0,.2,0,.2])
title('Diff by Active of Slope for FR=m*Ev+b')
xlabel('Differnce in Active Slopes: Low');
ylabel('Differnce in Active Slopes: Neutral');


set(gcf, 'Units', 'Inches')
set(gcf,'Position',[5.6167    3.6083    8.9667    3.6667]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
%saveas(gcf,['ParameterCompare.pdf'])  

%%
figure(9)
%histogram(MegaDataBase.LowDiff-MegaDataBase.NDiff)
hold on;
ofInt=(MegaDataBase_Neural.Pval_Center_Low>.05 & MegaDataBase_Neural.Pval_Center_Neutral>.05);
plot(MegaDataBase_Neural.LowDiff(ofInt),MegaDataBase_Neural.NDiff(ofInt),'ok');
NeitherNetralSig=sum(ofInt)

ofInt=(MegaDataBase_Neural.Pval_Center_Low>.05 & MegaDataBase_Neural.Pval_Center_Neutral<.05);
plot(MegaDataBase_Neural.LowDiff(ofInt),MegaDataBase_Neural.NDiff(ofInt),'or','MarkerFaceColor','r');
NeutralNetralSig=sum(ofInt)

ofInt=(MegaDataBase_Neural.Pval_Center_Low<.05 & MegaDataBase_Neural.Pval_Center_Neutral>.05);
plot(MegaDataBase_Neural.LowDiff(ofInt),MegaDataBase_Neural.NDiff(ofInt),'ob','MarkerFaceColor','b');
LowNetralSig=sum(ofInt)

ofInt=(MegaDataBase_Neural.Pval_Center_Low<.05 & MegaDataBase_Neural.Pval_Center_Neutral<.05);
plot(MegaDataBase_Neural.LowDiff(ofInt),MegaDataBase_Neural.NDiff(ofInt),'o','color',[0.6941    .5    0.8510],'MarkerFaceColor',[0.6941    .5    0.8510]);
BothNetralSig=sum(ofInt)

xlabel('Difference during Low (sp/sec)')
ylabel('Difference during Neutral (sp/sec)')
plot([0,9],[0,9],'--k');

legend({'Neither significant', 'Neutral significant','Low significant', 'Both significant'})

set(gcf, 'Units', 'Inches')
 set(gcf,'Position',[8.7000    1.0750    5.3000    4.3750]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,['NeutralCueComp_Final.pdf'])  
% %% Compare shift and H-rate effect
% 
% figure(8)
% subplot(1,2,1)
% hold on
% if sum((MegaDataBase.DiffInt_Low_Mem_p>.05))>0
% plot(MegaDataBase.Bias_Shift_POS(MegaDataBase.DiffInt_Low_Mem_p>.05), abs(MegaDataBase.DiffInt_Low_Mem(MegaDataBase.DiffInt_Low_Mem_p>.05)),'o','MarkerEdgeColor','k',...
%               'MarkerFaceColor','w',...
%               'LineWidth',1.5, 'DisplayName', 'Non-Sig')
% end
% if sum((MegaDataBase.DiffInt_Low_Mem_p<.05))>0
% plot(MegaDataBase.Bias_Shift_POS(MegaDataBase.DiffInt_Low_Mem_p<.05),abs(MegaDataBase.DiffInt_Low_Mem(MegaDataBase.DiffInt_Low_Mem_p<.05)),'o','MarkerEdgeColor','b',...
%               'MarkerFaceColor','b',...
%               'LineWidth',1.5, 'DisplayName', 'Low: p<.05')
% end
% legend('Autoupdate','off','Location', 'best')
% title('Diff by Active of Intercept for FR=m*Ev+b')
% xlabel( 'Bias Shift (Pos is optimal)');
% ylabel({'Abs Difference in Intercepts by H rate;','May not be in correct direction'});
% plot([0,0],[0,30],'--k')
% axis([0,3,0,30]);
% 
% subplot(1,2,2)
% hold on
% if sum((MegaDataBase.DiffSlope_Low_Mem_p>.05))>0
% plot(MegaDataBase.Bias_Shift_POS(MegaDataBase.DiffSlope_Low_Mem_p>.05), abs(MegaDataBase.DiffSlope_Low_Mem(MegaDataBase.DiffSlope_Low_Mem_p>.05)),'o','MarkerEdgeColor','k',...
%               'MarkerFaceColor','w',...
%               'LineWidth',1.5, 'DisplayName', 'Low: p>.05')
% end
% if sum((MegaDataBase.DiffSlope_Low_Mem_p<.05))>0
% plot(MegaDataBase.Bias_Shift_POS(MegaDataBase.DiffSlope_Low_Mem_p<.05),abs(MegaDataBase.DiffSlope_Low_Mem(MegaDataBase.DiffSlope_Low_Mem_p<.05)),'o','MarkerEdgeColor','b',...
%               'MarkerFaceColor','b',...
%               'LineWidth',1.5, 'DisplayName', 'Low: p<.05')
% end
% %legend('Autoupdate','off','Location', 'best')
% title('Diff by Active of Slope for FR=m*Ev+b')
% xlabel( 'Bias Shift (Pos is optimal)');
% ylabel({'Abs Difference in Slopes by H rate;','May not be in correct direction'});
% plot([0,0],[0,.2],'--k')
% axis([0,3,0,.2]);
%%

%
figure(10)
subplot(1,3,1)
hold on
if sum((MegaDataBase_Neural.T1_Diff_Mem>.05))>0
plot(1000*MegaDataBase_Neural.T1_GSFR_M(MegaDataBase_Neural.T1_Diff_Mem>.05 & MegaDataBase_Neural.T1_GSFR_M>0 & MegaDataBase_Neural.T1_NeutralFR_M>0), 1000*(MegaDataBase_Neural.T1_NeutralFR_M(MegaDataBase_Neural.T1_Diff_Mem>.05 & MegaDataBase_Neural.T1_GSFR_M>0 & MegaDataBase_Neural.T1_NeutralFR_M>0)),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Low: p>.05')
end

if sum((MegaDataBase_Neural.T1_Diff_Mem<.05))>0
plot(1000*MegaDataBase_Neural.T1_GSFR_M(MegaDataBase_Neural.T1_Diff_Mem<.05 & (MegaDataBase_Neural.T1_GSFR_M*1000)>0.5 & MegaDataBase_Neural.T1_NeutralFR_M>0),1000*(MegaDataBase_Neural.T1_NeutralFR_M(MegaDataBase_Neural.T1_Diff_Mem<.05 & (MegaDataBase_Neural.T1_GSFR_M*1000)>0.5 & MegaDataBase_Neural.T1_NeutralFR_M>0)),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5, 'DisplayName', 'Low: p<.05')
end
title('T1')
xlabel('Guided Saccade FR (sp/sec)')
ylabel('Neutral Hazard FR (sp/sec)')
plot([0,25],[0,25],'--k');
T1SigDiff=signrank(1000*MegaDataBase_Neural.T1_GSFR_M(MegaDataBase_Neural.T1_GSFR_M>0 & MegaDataBase_Neural.T1_NeutralFR_M>0)-1000*MegaDataBase_Neural.T1_NeutralFR_M(MegaDataBase_Neural.T1_GSFR_M>0 & MegaDataBase_Neural.T1_NeutralFR_M>0))

subplot(1,3,2)
hold on
if sum((MegaDataBase_Neural.N_Diff_Mem>.05))>0
plot(1000*MegaDataBase_Neural.N_GSFR_M(MegaDataBase_Neural.N_Diff_Mem>.05 & (1000*MegaDataBase_Neural.N_GSFR_M)>.5 & (1000*MegaDataBase_Neural.N_NeutralFR_M)>.5), 1000*(MegaDataBase_Neural.N_NeutralFR_M(MegaDataBase_Neural.N_Diff_Mem>.05 & (1000*MegaDataBase_Neural.N_GSFR_M)>.5 & (1000*MegaDataBase_Neural.N_NeutralFR_M)>.5)),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Low: p>.05')
end
if sum((MegaDataBase_Neural.N_Diff_Mem<.05))>0
plot(1000*MegaDataBase_Neural.N_GSFR_M(MegaDataBase_Neural.N_Diff_Mem<.05 & (1000*MegaDataBase_Neural.N_GSFR_M)>.5 & (1000*MegaDataBase_Neural.N_NeutralFR_M)>.5),1000*(MegaDataBase_Neural.N_NeutralFR_M(MegaDataBase_Neural.N_Diff_Mem<.05 & (1000*MegaDataBase_Neural.N_GSFR_M)>.5 & (1000*MegaDataBase_Neural.N_NeutralFR_M)>.5)),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5, 'DisplayName', 'Low: p<.05')
end
title('N')
xlabel('Guided Saccade FR (sp/sec)')
plot([0,25],[0,25],'--k');
NSigDiff=signrank(1000*MegaDataBase_Neural.N_GSFR_M(MegaDataBase_Neural.N_GSFR_M>0 & MegaDataBase_Neural.N_NeutralFR_M>0)-1000*MegaDataBase_Neural.N_NeutralFR_M(MegaDataBase_Neural.N_GSFR_M>0 & MegaDataBase_Neural.N_NeutralFR_M>0))

subplot(1,3,3)
hold on
if sum((MegaDataBase_Neural.T2_Diff_Mem>.05))>0
plot(1000*MegaDataBase_Neural.T2_GSFR_M(MegaDataBase_Neural.T2_Diff_Mem>.05 & (1000*MegaDataBase_Neural.T2_GSFR_M)>.5 & (1000*MegaDataBase_Neural.T2_NeutralFR_M)>.5), 1000*(MegaDataBase_Neural.T2_NeutralFR_M(MegaDataBase_Neural.T2_Diff_Mem>.05 & (1000*MegaDataBase_Neural.T2_GSFR_M)>.5 & (1000*MegaDataBase_Neural.T2_NeutralFR_M)>.5)),'o','MarkerEdgeColor','k',...
              'MarkerFaceColor','w',...
              'LineWidth',1.5, 'DisplayName', 'Low: p>.05')
end
if sum((MegaDataBase_Neural.T2_Diff_Mem<.05))>0
plot(1000*MegaDataBase_Neural.T2_GSFR_M(MegaDataBase_Neural.T2_Diff_Mem<.05 & (1000*MegaDataBase_Neural.T2_GSFR_M)>.5 & (1000*MegaDataBase_Neural.T2_NeutralFR_M)>.5),1000*(MegaDataBase_Neural.T2_NeutralFR_M(MegaDataBase_Neural.T2_Diff_Mem<.05 & (1000*MegaDataBase_Neural.T2_GSFR_M)>.5 & (1000*MegaDataBase_Neural.T2_NeutralFR_M)>.5 )),'o','MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5, 'DisplayName', 'Low: p<.05')
end
title('T2')
xlabel('Guided Saccade FR (sp/sec)')
T2SigDiff=signrank(1000*MegaDataBase_Neural.T2_GSFR_M(MegaDataBase_Neural.T2_GSFR_M>0 & MegaDataBase_Neural.T2_NeutralFR_M>0)-1000*MegaDataBase_Neural.T2_NeutralFR_M(MegaDataBase_Neural.T2_GSFR_M>0 & MegaDataBase_Neural.T2_NeutralFR_M>0))
plot([0,25],[0,25],'--k');


set(gcf, 'Units', 'Inches')
set(gcf,'Position',[0.3583    3.5521   15.6167    4.3750]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,['T1_N_T2_RdespCompare.pdf'])  


% close all
% figure
% hold on
% subplot(1,2,1);
% hold on;
% plot(lowH_Haz,'-or')
% plot(highH_Haz,'-ob');
% plot([0,length(lowH_Haz)],[.5,.5],'b')
% plot([0,length(lowH_Haz)],[.05,.05],'r')
% ylabel('Fit H')
% xlabel('Session (or Half Session)')
% set(gca, 'Fontsize', 14);
% subplot(1,2,2);
% hold on;
% plot(lowH_Log,'-or');
% plot(highH_Log,'-ob');
% set(gca, 'Fontsize', 14);
% 
% ylabel('Psychometric Curve Bias Term')
% xlabel('Session (or Half Session)')
%  set(gcf, 'Units', 'Inches')
%  
%      set(gcf,'Position',[ 1.0167    3.5333   12.8750    4.3750]);
% pos=get(gcf, 'Position');
% set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
% saveas(gcf,['Ci_Behavior Summary.pdf'])  

% %% Look for directionality 
% listOfSig=[];
% listOfNoMem=[];
% Index=[];
% Selectivity=[];
% totalneurons=0;
% for sessions=1:length(FitDirs)
% for neurons=1:size(FitDirs{1,sessions},1)
%         if sum(isnan(FitDirs{1,sessions}(neurons,:)))~=4
% for Periods=1:4
% if FitDirs{1,sessions}(neurons,Periods)<.05
% disp(FitDirs{2,sessions}{1,neurons}{1,1}{Periods,2})
% Selectivity=[Selectivity;FitDirs{2,sessions}{1,neurons}{1,1}{Periods,2}];
% Index=[Index;FitDirs{3,sessions}(neurons,Periods)];
% listOfSig=[listOfSig;sessions,neurons,Periods];
% 
% end
% end
% totalneurons=totalneurons+1
% 
%     else
%         listOfNoMem=[listOfNoMem;sessions,neurons];
% end
% end
% end
% % %% Look for H rate differences
% % HighHDiff=PeriodMeanFR{1,currentFile-1}{1,2}(4:6,2,1,1,1)-PeriodMeanFR{1,currentFile-1}{1,2}(4:6,2,2,1,1)
% % LowHDiff=PeriodMeanFR{1,currentFile-1}{1,2}(4:6,1,1,1,1)-PeriodMeanFR{1,currentFile-1}{1,2}(4:6,1,2,1,1)
% % mean(abs(HighHDiff))-mean(abs(LowHDiff))
