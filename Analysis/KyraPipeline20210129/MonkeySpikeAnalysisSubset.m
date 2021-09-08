function [AvgFRTotal,STDFRTotal, NumTri ] = MonkeySpikeAnalysisSubset(data, subby, figgy,b, TACP_Cuttoff)
%This function creates PSTH of data for AODR (and option of Raster)
%task separated out by:
%H=Hazard Rates
%A=Which target is Active
%S=Which subset you are interested in
%1= Correct vs Error Trials
%2= Monkey stayed with choice from what was correct last time or switched
%3= Monkey Chose T1 or T2




%%
datahold=data;
clear Mfr
clear Sfr
col={'r','--r','b','--b';'m','--m','c','--c'};

%First, need to split up by evidence,
% Put sample angles on -180<x<180, useful for when targets are 45/315,
% maybe not so good for 135/225
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
%Collect all the unique locations (the rounding is because I found the
%divisions made the samples not always line up whether T1 or T2 was
%accurate
Evidences=unique(round(datahold.Sample_Angle));
datahold.Sample_Angle=round(datahold.Sample_Angle);
%Remove any Nan's that may have snuck in
Evidences=Evidences(~isnan(Evidences));
Active=unique(datahold.Choice);
Active=Active(~isnan(Active));

NumTri=[];
clear fr
%If you are giving a file that is an aggregate and want to do by session
num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});
allFR=cell(1,2);
for neuro=1:unique(datahold.Num_Neuron);
    fig=figure(figgy+neuro-1);
    tiledlayout(3,3, 'TileSpacing', 'compact')
    spiketimes=datahold{:,24+neuro};
    %Predefine some varaibles
    for Ev=1:length(Evidences)
        for H=1:2
            fr{Ev,H}=[];
            Rastas{Ev}=[];
            TotalBinnedSpikes{Ev,H}=[];
        end
    end
    
    Hs=unique(datahold.H_Rate);
    
    %Define what subset you are interested in
    if subby==1 %corr vs err
        opts=[1,0];
        catagory=datahold.Correct;
        namer={'Corr','Err'};
        label='Corr';
    elseif subby==2 %sw vs not
        opts=[0,1];
        catagory=datahold.Switch;
        namer={'Stay', 'Switch'};
        label='Switch';
    elseif subby==3  %L vs R
        opts=Active;
        catagory=datahold.Choice;
        namer={'Choose1','Choose2'};
        label='Choice';
    end
    
    for Sess=1%:num_sessions
        %Find out which entries in the mega file are from a particular
        %session
        sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
        for H=1:length(Hs)  %For each H rate
            for A=1:2 %Which is the ACTIVE Target
                for S=1:2  % Which subset factor (e.g. corr/err)
                    for Ev=1:length(Evidences)
                        TotalBinnedSpikes{Ev}=[];
                        
                        if H==1
                            % Specifically interested in if there is a
                            % difference between T1 active and T2 active
                            % for Low H (H=1), so excluding trials around
                            % the switch
                            subsettle= datahold.Active_Angle==Active(A) & datahold.TACP_Actual>TACP_Cuttoff;
                        else
                            subsettle= datahold.Active_Angle==Active(A);
                            
                        end
                        
                        % Pick the trials from thes session, that are of
                        % the correct sample angle (Ev), that are of the
                        % correct H rate, that are of the correct active
                        % target, and are of the correct subset catagory
                        % (err/corr eg)
                        TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1 & catagory==opts(S));
                        NumTri(Ev,(4*(H-1)+2*(A-1)+S))=length(TOI);
                        
                        %From 500 ms before stim offset to 3 seconds after
                        %in bins of 100 ms
                        binsize=100;
                        bins=[-500:binsize:3000];
                        for trials=1:length(TOI)
                            
                            %Raster
                            spiky=spiketimes{TOI(trials)};
                            spiky=spiky(spiky<4000)-datahold.SampleOff(TOI(trials));
                            % If took less than 5 sec to fixate, and made
                            % it to the end of trial (i.e. FP went off)
                            wantRasta=0;
                            if datahold.SampleOn(TOI(trials))<5000 & datahold.SampleOff(TOI(trials))<5000 & datahold.FPOff(TOI(trials)) & wantRasta
                                %Need to organize the figures here, but
                                %it's a lot (2H*2A*2S*9Ev=72 Rasters!)
                                hold on
                                scatter(spiky,trials*ones(1,length(spiky)),'.k');
                                scatter(datahold.SampleOn(TOI(trials)),trials,'.g');
                                scatter(datahold.FPOff(TOI(trials)),trials,'.r');
                            end
                            
                            
                            %PSTH
                            binnedspikes=[];
                            for bin=1:length(bins)-2
                                binnedspikes(1,bin)=sum(spiky>bins(bin) & spiky<bins(bin+2));
                            end
                            
                            %All the number of spikes in bin across Ev an
                            TotalBinnedSpikes{Ev}=[TotalBinnedSpikes{Ev};binnedspikes];
                            FPOfftime=datahold.FPOff(TOI(trials))-datahold.SampleOff(TOI(trials));
                        end
                        
                        
                        if ~isempty(TotalBinnedSpikes{Ev})
                            AvgFRPart=nanmean(TotalBinnedSpikes{Ev},1)/binsize*1000;
                            STDFRPart=nanstd(TotalBinnedSpikes{Ev},1)/binsize*1000;
                        else
                            AvgFRPart=nans(1,length(bins)-2);
                            STDFRPart=nans(1,length(bins)-2);;
                        end
                        
                        %Start to plot Avg FR
                        
                        hold on
                        nexttile(Ev)
                        
                        %Optional, can sutract baseline=mean FR of
                        %pre-Sample Off period
                        baseline=0; %nanmean(MNSpkCnt(1:4));
                        hold on
                        plot(bins(2:end-1),AvgFRPart-baseline, col{H,S+2*(A-1)});
                        
                        title({['Sample Angle= ', num2str(Evidences(Ev))];['Neuron ' datahold.Properties.VariableNames{24+neuro}]},'Interpreter', 'none');
                        
                        ylabel('Avg FR');
                        xlabel('Time(ms) From Samp Off');
                        
                        %Store AvgFR
                        AvgFRTotal{Ev,H,A,S,neuro}=AvgFRPart;
                        STDFRTotal{Ev,H,A,S,neuro}=STDFRPart;
                    end
                    
                end
            end
        end
        
        legend(['H= ', num2str(Hs(1)), ' T1 Active ',namer{1}], ['H= ', num2str(Hs(1)),' T1 Active ',namer{2}],['H= ', num2str(Hs(1)),' T2 Active ',namer{1}], ['H= ', num2str(Hs(1)),' T2 Active ',namer{2}],['H= ', num2str(Hs(2)), ' T1 Active ',namer{1}], ['H= ', num2str(Hs(2)),' T1 Active ',namer{2}],['H= ', num2str(Hs(2)),' T2 Active ',namer{1}], ['H= ', num2str(Hs(2)),' T2 Active ',namer{2}],'location','eastoutside','AutoUpdate','off');
        
        %Put a marker of when FP off on all plots
        for Ev=1:9
            nexttile(Ev)
            plot([FPOfftime,FPOfftime],[0, 10],'--k');
        end
    end
 %Save the plot
    set(gcf, 'Units', 'Inches')
     set(gcf,'Position',[1.3083    0.6333   13.9167    7.5083]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,[b,'_neuron_', num2str(neuro), '_PETH_',label,'.pdf'])  
    if neuro==datahold.Num_Neuron
        NumTri=num2cell(NumTri);
        NumTri=cell2table(NumTri,'VariableNames',{['HL_A1_',namer{1}],['HL_A1_',namer{2}],['HL_A2_',namer{1}],['HL_A2_',namer{2}],['HH_A1_',namer{1}],['HH_A1_',namer{2}],['HH_A2_',namer{1}],['HH_A2_',namer{2}]},'Rownames',{['Ev_',num2str(Evidences(1))],['Ev_',num2str(Evidences(2))],['Ev_',num2str(Evidences(3))],['Ev_',num2str(Evidences(4))],['Ev_',num2str(Evidences(5))],['Ev_',num2str(Evidences(6))],['Ev_',num2str(Evidences(7))],['Ev_',num2str(Evidences(8))],['Ev_',num2str(Evidences(9))]});
    end
end
end
