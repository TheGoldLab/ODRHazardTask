%This script will run the behavioral analysis for aODR trials for all the
%.mat files in data_path.  The data will be plotted per session and put in
%the results_path. All data of the analysis will be saved in
%'Aggro_Behav.mat' in the results_path under the variable MegaDataBase.
        %MegaDataBase base contains the following varriables
            %1=File
            %2=Monkey
            %3=Data
            %4=Difference in Subjective Hazard rates
            %5=Difference in fit logistic bias terms
            %6=Lapse rate from the logistic fit
            %7=table of all the raw fits for both normative model and logistif fit
                %Row 1=low H rate
                %Row 2=neutral H rate
                %Columns: Fit means from the Normative model fit, Log means from the
                %logistic fit.

%The script will first lookd for a file called "Aggro_Behav.mat" in the
%results_path and load it.  It will then see how many sessions


% Where are the files you are running
data_path = 'C:\Users\kyras\Desktop\UPenn 17-18\Matlab_AODR\ODRHazardTask-main\BehOnly';%AVData/794'; % DatatoParse2';
results_path = 'C:\Users\kyras\Desktop\UPenn 17-18\Matlab_AODR\ODRHazardTask-main\BehOnly\Results';%AVData/794'; % DatatoParse2';
cd(data_path)
mkdir(results_path)

%Collect the .mat files in the directory
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);

%Set up the database
MegaDataBase=table;
temptoAdd=table;
close all
cd(results_path)
%%
if exist('Aggro_Behav.mat','file')==2
    
    load('Aggro_Behav.mat');
    newfiles=1;
    for currentFile =  1:nfiles
        
        cd(data_path);
        
        %(responses within this many degree of target are considered a response, rest are NaN)
        ResponseCutoff=30;
        theFile= files(currentFile).name;
        if sum(ismember(MegaDataBase.File,theFile))==0
            
            temp=load(theFile);
            data=temp.data;
            [~,b,~]=fileparts(theFile);
            cd(data_path);
            
            %Parse the files for aODR trials, will gather all the
            %relevant information
            %Response cut off: Will reject saccades to >Response Cutt off away from
            %a target
            %2= aODR trials (1 will do ODR trials)
            %data=data from the FIRA
            %theFile=name of the file
            TaskTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,2,data,theFile,data_path);
            
            %Start making figures/doing analysis
            
            %For low H rate, how many trials after the CP do you want to consider,
            %useful when looking at neural things, but not used for behavior;
            TACP_Cuttoff=0;
            if ~isempty(TaskTrials)
                [Fits,TrialNum]=FindHRateQuant(TaskTrials,1+newfiles,b,1,results_path);
                temptoAdd{newfiles,1}=string(theFile); %file name
                temptoAdd{newfiles,2}={theFile(1:2)}; %Monkey
                temptoAdd{newfiles,3}={theFile(4:13)}; %Date
                temptoAdd{newfiles,4}=Fits.Fit_H(2)-Fits.Fit_H(1); %Hrate_Shift
                temptoAdd{newfiles,5}=Fits.LogBias(1)-Fits.LogBias(2);
                temptoAdd{newfiles,6}=Fits.LogLaps(1);
                temptoAdd{newfiles,7}={Fits};
                newfiles=newfiles+1;
            end
        end
    end
    if newfiles>1
        temptoAdd.Properties.VariableNames={'File','Monkey','Date','Delta_H','Delta_Bias','Lapse_Log','All_FitData'};

    MegaDataBase=[MegaDataBase;temptoAdd];
    MegaDataBase = sortrows(MegaDataBase,'File','ascend');
    cd (results_path)
    save('Aggro_Behav.mat','MegaDataBase');
    end
else
    
    for currentFile =  1:nfiles
        
        cd(data_path);
        
        %(responses within this many degree of target are considered a response, rest are NaN)
        ResponseCutoff=30;
        theFile= files(currentFile).name;
        temp=load(theFile);
        data=temp.data;
        [~,b,~]=fileparts(theFile);
        cd(data_path);
        
        %Parse the files for aODR trials, will gather all the
        %relevant information
        %Response cut off: Will reject saccades to >Response Cutt off away from
        %a target
        %2= aODR trials (1 will do ODR trials)
        %data=data from the FIRA
        %theFile=name of the file
        TaskTrials=KAS_Parser_794_DelayNeural(ResponseCutoff,2,data,theFile,data_path);
        
        %Start making figures/doing analysis
        
        %For low H rate, how many trials after the CP do you want to consider,
        %useful when looking at neural things, but not used for behavior;
        TACP_Cuttoff=0;
        if ~isempty(TaskTrials)
            [Fits,TrialNum]=FindHRateQuant(TaskTrials,1+currentFile,b,1,results_path);
            MegaDataBase{currentFile,1}=string(theFile); %file name
            MegaDataBase{currentFile,2}={theFile(1:2)}; %Monkey
            MegaDataBase{currentFile,3}={theFile(4:13)}; %Date
            MegaDataBase{currentFile,4}=Fits.Fit_H(2)-Fits.Fit_H(1); %Hrate_Shift
            MegaDataBase{currentFile,5}=Fits.LogBias(1)-Fits.LogBias(2);
            MegaDataBase{currentFile,6}=Fits.LogLaps(1);
            MegaDataBase{currentFile,7}={Fits};
        end
    end
    MegaDataBase.Properties.VariableNames={'File','Monkey','Date','Delta_H','Delta_Bias','Lapse_Log','All_FitData'};
    
    cd (results_path)
    save('Aggro_Behav.mat','MegaDataBase');
end


%%
monklist=unique(MegaDataBase.Monkey);
%SessionOfInterest
SOI=[];
for monk=1:length(unique(MegaDataBase.Monkey))
 
    
    for i=1:height(MegaDataBase)
        SOI(i,1)=strcmp(MegaDataBase.Monkey(i),monklist(monk));
    end
    stdBias(monk)=nanstd(MegaDataBase.Delta_Bias(find(SOI)));
    stdH(monk)=nanstd(MegaDataBase.Delta_H(find(SOI)));
    working=MegaDataBase(find(SOI),:)
    % For if you want to include sessions with values >3 STD away from mean
    % & abs(MegaDataBase.Delta_Bias)<(3*(stdBias(monk))) & abs(MegaDataBase.Delta_H)<(3*stdH(monk)),:);
    
    [a,b]=unique(working.Delta_Bias);
    
    sortedb=sort(b(~isnan(a)));
    
     figure(1)
     %Plot the difference in Subjective H rate
    subplot(3,length(monklist),monk+2*(length(monklist)))
    plot(working.Delta_H(sortedb),'-ok','Linewidth',2)
    hold on
    plot([0,length(working.Delta_H(sortedb))],[0,0],'--k')
    if monk==1
        ylabel({'{\Delta}Subjective H rate'})
    end
    title({'H rate Difference over time'})
    
    %Plot difference in fit Bias term
    subplot(3,length(monklist),monk+(length(monklist)))
    plot(working.Delta_Bias(sortedb),'-ok','Linewidth',2)
    hold on
    plot([0,length(working.Delta_Bias(sortedb))],[0,0],'--k')
    if monk==1
        ylabel({'{\Delta}Bias'})
    end
    
    title('Bias Difference over time')
    set(gcf, 'Units', 'Inches')
    %set(gcf,'Position',[7.8000    1.9167    6.8583    6.2250]);
    pos=get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
   
    %Some statistics
    sigdiffH(monk)=signrank(working.Delta_H(sortedb))
    sigdiffBias(monk)=signrank(working.Delta_Bias(sortedb))
    avgdiffH(monk)=mean(working.Delta_H(sortedb));
    avgdiffBias(monk)=mean(working.Delta_Bias(sortedb))
    
    %Plot Lapse rate over sessions
    subplot(3,length(monklist),monk)
    plot(working.Lapse_Log(sortedb),'-ok','Linewidth',2)
    hold on
    plot([0,length(working.Delta_H(sortedb))],[0,0],'--k')
    title({monklist{monk};'Lapse Rate'});
    xlabel({'Session Number'})
    if monk==1
        ylabel('Lapse Rate')
    end
    LapseMean(monk)=nanmean(working.Lapse_Log(sortedb));
end
saveas(gcf,['BehavioralMeasures_2022.pdf'])


%%