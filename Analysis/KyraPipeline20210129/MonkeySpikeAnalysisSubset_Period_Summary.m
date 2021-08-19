function [ Mfr,Sfr ] = MonkeySpikeAnalysisSubset_Period_Summary(data, subby, figgy,b,TACP_Cutoff,Alice)
%This function Avg Fr as fx of 
%H=Hazard Rates
%A=Which target is Active
%S=Which subset you are interested in
    %1= Correct vs Error Trials
    %2= Monkey stayed with choice from what was correct last time or switched
    %3= Monkey Chose T1 or T2
%For each trial Period:
    %1= Visual
    %2= Memory
    %3= Saccade
    %4= Reward period?


datahold=data;
col={'-or','--dr','-ob','--db';'-om','--dm','-oc','--dc'};
colsBar={[1,0,0],[1,.5,.5],[0,0,1], [.5,.5,1],[1,0,0],[0,0,1],[1,0,1], [0,1,1]};


%First, need to split up by evidence, maybe later also by if chose A or B
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
Evidences=unique(round(datahold.Sample_Angle));
datahold.Sample_Angle=round(datahold.Sample_Angle);
Evidences=Evidences(~isnan(Evidences));
Choices=unique(datahold.Choice);
Choices=Choices(~isnan(Choices));

NumTri=[];
clear fr
num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});
allFR=cell(1,2);
for neuro=1:unique(datahold.Num_Neuron);
    spiketimes=datahold{:,24+neuro};
    for Ev=1:length(Evidences)
        for H=1:2
            Mfr=[];
            Sfr=[];
            Baselist=[];
        end
    end
    
    Hs=[.05,.50];
    
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
        opts=Choices;
        catagory=datahold.Choice;
        namer={'Choose1','Choose2'};
        label='Choice';
    end
    
    
                    
                    figr=figure(figgy+neuro-1);
                    if ~Alice
    tiledlayout(1,5, 'TileSpacing', 'compact')    ;
                    else
                          tiledlayout(5,8, 'TileSpacing', 'compact')    ;
                    end
    Sess=1;%:num_sessions
    %Find out which entries in the mega file are from a particular
    %session
    sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
    for H=1:2
        for A=1:2 %Which is the Active target
            for S=1:2  % which of subset options
                for P=1:5
                    if P==1
                        starter=datahold.SampleOn;%datahold.SacOn; %datahold.SampleOff;
                        ender=datahold.SampleOff; %starter+500;%datahold.FPOff;
                        titly={[fileNAMES{Sess}],['Visual'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]]};
                    elseif P==2
                        starter=datahold.SampleOff;%datahold.SacOn; %datahold.SampleOff;
                        ender=datahold.FPOff; %starter+500;%datahold.FPOff;
                        titly={[fileNAMES{Sess}],['Memory'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]]};
                    elseif P==3
                        starter=datahold.FPOff; %datahold.SampleOff;
                        ender=datahold.TargAcq;%datahold.FPOff;
                        titly={[fileNAMES{Sess}],['Motor'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]];};
                   
                    elseif P==4
                        starter=datahold.TargAcq; %datahold.SampleOff;
                        ender=datahold.TargAcq+500;%datahold.FPOff;
                        titly={[fileNAMES{Sess}],['Reward'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]];};
                    else
                        starter=datahold.TargAcq+500; %datahold.SampleOff;
                        ender=datahold.TargAcq+1500;%datahold.FPOff;
                        titly={[fileNAMES{Sess}],['PostRew'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]];};
                        
                    end

                    for Ev=1:length(Evidences)
                         fr{Ev}=[];

                        if H==1
                            subsettle= datahold.Active_Angle==Choices(A) & datahold.TACP_Actual>TACP_Cutoff; % Active_Angle
                        else
                            subsettle= datahold.Active_Angle==Choices(A);
                        end
                        TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1 & catagory==opts(S));
                        NumTri(Ev,(4*(H-1)+2*(A-1)+S))=length(TOI);
                        if length(TOI)>1
                            for trials=1:length(TOI)
                                spks=sum([spiketimes{TOI(trials)}]>starter(TOI(trials)) & [spiketimes{TOI(trials)}]<ender(TOI(trials)));
                                fr{Ev}=[fr{Ev};spks/(ender(TOI(trials))-starter(TOI(trials)))];
                                Baselist=[Baselist,sum([spiketimes{TOI(trials)}]<[datahold.SampleOn(TOI(trials))])/[datahold.SampleOn(TOI(trials))]];
                            end
                            
                            Mfr(Ev,H,A,S,P)=nanmean(fr{Ev}*1000);
                            Sfr(Ev,H,A,S,P)=nanstd(fr{Ev}*1000);
                        else
                            Mfr(Ev,H,A,S,P)=NaN;
                            Sfr(Ev,H,A,S,P)=NaN;
                        end
                    end

                   

                    if~ Alice
                         nexttile(P);
                        plot(Evidences,Mfr(:,H,A,S,P),col{H,S+2*(A-1)},'Markersize',5);
                                        ylabel('FR (spikes/sec)')
                    xlabel('Sample Angle (deg)')
                    hold on
                    title(titly,'Interpreter', 'none');
                    else
                        nexttile(8*(P-1)+S+2*(A-1)+4*(H-1))
                        bar(Evidences,Mfr(:,H,A,S,P));
                        if 8*(P-1)+S+2*(A-1)+4*(H-1)==1
                            ylabel('Visual FR')
                            title([titly, 'H=.05, A1, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==2
                            title([titly, 'H=.05, A1, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==3
                            title([titly, 'H=.05, A2, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==4
                            title([titly, 'H=.05, A2, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==5
                            title([titly, 'H=.5, A1, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==6
                            title([titly, 'H=.5, A1, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==7
                            title([titly, 'H=.5, A2, S=', num2str(namer{S})],'Interpreter','none');
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==8
                            title([titly, 'H=.5, A2, S=', num2str(namer{S})],'Interpreter','none');

                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==9
                            ylabel('Memory FR')
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==17
                            ylabel('Motor FR')
                        elseif 8*(P-1)+S+2*(A-1)+4*(H-1)==25
                            ylabel('Reward FR')
                        end
                        
                        
                    end

                end

            end
        end
    end
    
   
  if ~Alice
    legend(['H=.05, A1; ', namer{1}], ['H=.05, A1; ', namer{2}], ['H=.05, A2; ', namer{1}], ['H=.05, A2; ', namer{2}],['H=.5, A1; ', namer{1}], ['H=.5, A1; ', namer{2}], ['H=.5, A2; ', namer{1}], ['H=.5, A2; ', namer{2}], 'Location', 'eastoutside')
  end
    
    set(gcf, 'Units', 'Inches')
    if ~Alice
     set(gcf,'Position',[ 1.3083    4.4917   13.9167    3.6499]);
    else
         set(gcf,'Position',[1.3083    0.4250   13.1834    7.7166]);
         
    end
    
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
if ~Alice
saveas(gcf,[b,'_neuron_', num2str(neuro), '_Period_Summary_',label,'.pdf'])  
else
    saveas(gcf,[b,'_neuron_', num2str(neuro), '_Period_Summary_Bar_',label,'.pdf'])  
end
end

if neuro==datahold.Num_Neuron
    NumTri=num2cell(NumTri);
    NumTri=cell2table(NumTri,'VariableNames',{['HL_A1_',namer{1}],['HL_A1_',namer{2}],['HL_A2_',namer{1}],['HL_A2_',namer{2}],['HH_A1_',namer{1}],['HH_A1_',namer{2}],['HH_A2_',namer{1}],['HH_A2_',namer{2}]},'Rownames',{['Ev_',num2str(Evidences(1))],['Ev_',num2str(Evidences(2))],['Ev_',num2str(Evidences(3))],['Ev_',num2str(Evidences(4))],['Ev_',num2str(Evidences(5))],['Ev_',num2str(Evidences(6))],['Ev_',num2str(Evidences(7))],['Ev_',num2str(Evidences(8))],['Ev_',num2str(Evidences(9))]});
end
end

