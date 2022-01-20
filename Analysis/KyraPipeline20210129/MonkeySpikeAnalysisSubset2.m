function [AvgFRTotal,STDFRTotal, NumTri ] = MonkeySpikeAnalysisSubset2(data)
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
col={'r','--r','b','--b','g', 'g--',;'m','--m','c','--c','y', 'y--'};

%First, need to split up by Lidence,
% Put sample angles on -180<x<180, useful for when targets are 45/315,
% % maybe not so good for 135/225
% datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
% %Collect all the unique locations (the rounding is because I found the
% %divisions made the samples not always line up whether T1 or T2 was
% %accurate
% Lidences=unique(round(datahold.Sample_Angle));
% datahold.Sample_Angle=round(datahold.Sample_Angle);
% %Remove any Nan's that may have snuck in
% Lidences=Lidences(~isnan(Lidences));
% Active=unique(datahold.Choice);
% Active=Active(~isnan(Active));

LLR_Ch_tr = round(datahold.LLR_for_Change,3);
LLR_Ch = unique(LLR_Ch_tr(2:end)); %first trial is nan

NumTri=[];
clear fr
%If you are giving a file that is an aggregate and want to do by session
num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});
allFR=cell(1,2);

AvgFRTotal = cell(length(LLR_Ch),length(unique(datahold.H_Rate)),unique(datahold.Num_Neuron));
STDFRTotal = cell(length(LLR_Ch),length(unique(datahold.H_Rate)),unique(datahold.Num_Neuron));
for neuro=1:unique(datahold.Num_Neuron)
%     fig=figure(figgy+neuro-1);
    fig = figure;
    tiledlayout(3,3, 'TileSpacing', 'compact')
    spiketimes=datahold{:,24+neuro};
    %Predefine some varaibles
    for L = 1:length(LLR_Ch) % L=1:length(Lidences)
        for H=1:2
            fr{L,H}=[];
            Rastas{L}=[];
            TotalBinnedSpikes{L,H}=[];
        end
    end
    
    Hs=unique(datahold.H_Rate);
    
    %Define what subset you are interested in
%     if subby==1 %corr vs err
%         opts=[1,0];
%         catagory=datahold.Correct;
%         namer={'Corr','Err'};
%         label='Corr';
%     elseif subby==2 %sw vs not
%         opts=[0,1];
%         catagory=datahold.Switch;
%         namer={'Stay', 'Switch'};
%         label='Switch';
%     elseif subby==3  %L vs R
%         opts=Active;
%         catagory=datahold.Choice;
%         namer={'Choose1','Choose2'};
%         label='Choice';
%     end
    
    for Sess=1%:num_sessions
        %Find out which entries in the mega file are from a particular
        %session
        sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
        for H=1:length(Hs)  %For each H rate
%             for A=1:2 %Which is the ACTIVE Target
%                 for S=1:2  % Which subset factor (e.g. corr/err)
                    for L=1:length(LLR_Ch)
                        TotalBinnedSpikes{L}=[];
                        
                        
%                         if H==1
%                             % Specifically interested in if there is a
%                             % difference between T1 active and T2 active
%                             % for Low H (H=1), so excluding trials around
%                             % the switch
%                             subsettle= datahold.Active_Angle==Active(A) & datahold.TACP_Actual>TACP_Cuttoff;
%                         else
%                             subsettle= datahold.Active_Angle==Active(A);
%                             
%                         end
                        
                        % Pick the trials from thes session, that are of
                        % the correct sample angle (L), that are of the
                        % correct H rate, that are of the correct active
                        % target, and are of the correct subset catagory
                        % (err/corr eg)
                        TOI=find(sessionBreakdown & LLR_Ch_tr==LLR_Ch(L) & abs(datahold.H_Rate-Hs(H))<.0001);
                        NumTri(L,H)=length(TOI);
                        
                        %From 500 ms before stim offset to 3 seconds after
                        %in bins of 100 ms
                        binsize=100;
                        bins=[-500:binsize:3000];
                        for tr=1:length(TOI)
                            
                            %Raster
                            spiky=spiketimes{TOI(tr)};
                            spiky=spiky(spiky<4000)-datahold.SampleOn(TOI(tr));
                            % If took less than 5 sec to fixate, and made
                            % it to the end of trial (i.e. FP went off)
                            wantRasta=0;
                            if datahold.SampleOn(TOI(tr))<5000 & datahold.SampleOn(TOI(tr))<5000 & datahold.FPOff(TOI(tr)) & wantRasta
                                %Need to organize the figures here, but
                                %it's a lot (2H*2A*2S*9L=72 Rasters!)
                                hold on
                                scatter(spiky,tr*ones(1,length(spiky)),'.k');
                                scatter(datahold.SampleOn(TOI(tr)),tr,'.g');
                                scatter(datahold.FPOff(TOI(tr)),tr,'.r');
                            end
                            
                            
                            %PSTH
                            binnedspikes=[];
                            for bin=1:length(bins)-2
                                binnedspikes(1,bin)=sum(spiky>bins(bin) & spiky<bins(bin+2));
                            end
                            
                            %All the number of spikes in bin across L an
                            TotalBinnedSpikes{L}=[TotalBinnedSpikes{L};binnedspikes];
                            FPOfftime=datahold.FPOff(TOI(tr))-datahold.SampleOn(TOI(tr));
                        end
                        
                        
                        if ~isempty(TotalBinnedSpikes{L})
                            AvgFRPart=nanmean(TotalBinnedSpikes{L},1)/binsize*1000;
                            STDFRPart=nanstd(TotalBinnedSpikes{L},1)/binsize*1000;
                        else
                            AvgFRPart=nans(1,length(bins)-2);
                            STDFRPart=nans(1,length(bins)-2);;
                        end
                        
                        %Start to plot Avg FR
                        
                        hold on
                        nexttile(L)
                        
                        %Optional, can sutract baseline=mean FR of
                        %pre-Sample Off period
                        baseline=0; %nanmean(MNSpkCnt(1:4));
                        hold on
                        plot(bins(2:end-1),AvgFRPart-baseline, col{2,H+H-1},'DisplayName',['H= ', num2str(Hs(H))]);
                        
                        title({['LLR_Change= ', num2str(LLR_Ch(L))];['Neuron ' datahold.Properties.VariableNames{24+neuro}(end-3:end)]},'Interpreter', 'none');
                        
                        ylabel('Avg FR');
                        xlabel('Time(ms) From Samp On');
                        
                        %Store AvgFR
%                         AvgFRTotal{L,H,A,S,neuro}=AvgFRPart;
%                         STDFRTotal{L,H,A,S,neuro}=STDFRPart;
                        AvgFRTotal{L,H,neuro}=AvgFRPart;
                        STDFRTotal{L,H,neuro}=STDFRPart;
                       
                    end
%                     
%                 end
%             end
        end
        
%         legend(['H= ', num2str(Hs(1)), ' T1 Active ',namer{1}], ['H= ', num2str(Hs(1)),' T1 Active ',namer{2}],['H= ', num2str(Hs(1)),' T2 Active ',namer{1}], ['H= ', num2str(Hs(1)),' T2 Active ',namer{2}],['H= ', num2str(Hs(2)), ' T1 Active ',namer{1}], ['H= ', num2str(Hs(2)),' T1 Active ',namer{2}],['H= ', num2str(Hs(2)),' T2 Active ',namer{1}], ['H= ', num2str(Hs(2)),' T2 Active ',namer{2}],'location','eastoutside','AutoUpdate','off');
%         legend(['H= ', num2str(Hs(1))], ['H= ', num2str(Hs(2))] );
        legend TOGGLE ON
        %Put a marker of when FP off on all plots
        for L=1%:9
%             nexttile(L)
            plot([FPOfftime,FPOfftime],[0, 10],'--k');
        end
    end
 %Save the plot
    set(gcf, 'Units', 'Inches')
     set(gcf,'Position',[1.3083    0.6333   13.9167    7.5083]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,[datahold.Session{1,:}(1:13),'_neuron_', num2str(neuro), '_PETH_','LLR_Change','.pdf'])  
%     if neuro==datahold.Num_Neuron
%         NumTri=num2cell(NumTri);
%         NumTri=cell2table(NumTri,'VariableNames',{['HL_A1_',namer{1}],['HL_A1_',namer{2}],['HL_A2_',namer{1}],['HL_A2_',namer{2}],['HH_A1_',namer{1}],['HH_A1_',namer{2}],['HH_A2_',namer{1}],['HH_A2_',namer{2}]},'Rownames',{['L_',num2str(Lidences(1))],['L_',num2str(Lidences(2))],['L_',num2str(Lidences(3))],['L_',num2str(Lidences(4))],['L_',num2str(Lidences(5))],['L_',num2str(Lidences(6))],['L_',num2str(Lidences(7))],['L_',num2str(Lidences(8))],['L_',num2str(Lidences(9))]});
%     end

% for neuro = 1:unique(datahold.Num_Neuron)
    figure
    hold on
    for H=1:length(Hs) 
        NormFR = nan(length(LLR_Ch),length(AvgFRPart));
        FRMat = vertcat(AvgFRTotal{:,H,neuro});
        FRMat(FRMat==0)=0.0000000001;
        NormFR = (FRMat-repmat(nanmean(FRMat),length(LLR_Ch),1))./repmat(nanmean(FRMat),length(LLR_Ch),1); %(std(FRMat)./sqrt(length(LLR_Ch)));

        plot(LLR_Ch,nanmean(NormFR(:,1:find(bins<FPOfftime,1,'last')),2),col{2,H+H-1},'DisplayName',['H= ', num2str(Hs(H))],'LineWidth',2)
        ylim([-1,1])
        xlim([-4,4])
        title({'LLR Change vs Normalized FR';fileNAMES{1}(1:13);['Neuron ' datahold.Properties.VariableNames{24+neuro}(end-3:end)]},'Interpreter', 'none');

        ylabel('% Diff from Avg LLR Trial');
        xlabel('LLR for Change');

    end
    legend TOGGLE ON
    set(gcf, 'Units', 'Inches')
    set(gcf,'Position',[7.8000    1.9167    6.8583    6.2250]);
    pos=get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
    saveas(gcf,[datahold.Session{1,:}(1:13),'Neuron ' datahold.Properties.VariableNames{24+neuro}(end-3:end),'_NeurometricFxn.pdf']) 
% end

end
end
