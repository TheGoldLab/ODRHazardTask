function [ Mfr,Sfr,HrateEffect, Comparison,SignRankLowInt,SignRankLowSlope,SignRankNInt,SignRankNSlope,IntRealLowDiff,IntRealNDiff,SlpeRealLowDiff,SlpeRealNDiff,pval ] = MonkeySpikeAnalysisSubset_Period_Summary(data, subby, figgy,b,TACP_Cutoff,Alice)
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
col={'or','--dr','ob','--db';'om','--dm','oc','--dc'};
col2={'r','--dr','b','--db';'m','--dm','c','--dc'};

colsBar={[1,0,0],[1,.5,.5],[0,0,1], [.5,.5,1],[1,0,0],[0,0,1],[1,0,1], [0,1,1]};


%First, need to split up by evidence, maybe later also by if chose A or B
if ~strncmp('Ci',data.Session{1},2)
    datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
end
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
    y=[];
    X=[];
    for Ev=1:length(Evidences)
        for H=1:2
            Mfr{neuro}=[];
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
                        fr{Ev,A}=[];
                          fr2{Ev,A}=[];
                        if H==1
                            subsettle= datahold.Active_Angle==Choices(A) & datahold.TACP_Actual>TACP_Cutoff; % Active_Angle
                        else
                            subsettle= datahold.Active_Angle==Choices(A);
                        end
                        TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1);% & catagory==opts(S));
                        %Want to have all trials for Neutral H rate neutral
                        %cue
                        TOI2=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1);
                        NumTri(Ev,(4*(H-1)+2*(A-1)+S))=length(TOI);
                        if length(TOI)>1
                            for trials=1:length(TOI)
                                spks=sum([spiketimes{TOI(trials)}]>starter(TOI(trials)) & [spiketimes{TOI(trials)}]<ender(TOI(trials)));
                                fr{Ev,A}=[fr{Ev,A};spks/(ender(TOI(trials))-starter(TOI(trials)))];

                                Baselist=[Baselist,sum([spiketimes{TOI(trials)}]<[datahold.SampleOn(TOI(trials))])/[datahold.SampleOn(TOI(trials))]];
                                if ~isnan(spks/(ender(TOI(trials))-starter(TOI(trials)))) ; %&S==1
                                    y=[y;spks/(ender(TOI(trials))-starter(TOI(trials)))*1000];
                                    X=[X;[Evidences(Ev),H,A,P]];
                                end
                                for trials=1:length(TOI2)
                                      spks2=sum([spiketimes{TOI2(trials)}]>starter(TOI2(trials)) & [spiketimes{TOI2(trials)}]<ender(TOI2(trials)));
                                fr2{Ev,A}=[fr2{Ev,A};spks2/(ender(TOI2(trials))-starter(TOI2(trials)))];

                                end
                                
                            end
                            
                            %Grab comparison info
                          
                            Mfr{neuro}(Ev,H,A,S,P)=nanmean(fr{Ev,A}*1000);
                            
                            Sfr(Ev,H,A,S,P)=nanstd(fr{Ev,A}*1000);
                            if Ev==5 && H==2
                                Mfr{neuro}(Ev,H,A,S,P)=nanmean(fr2{Ev,A}*1000);
                            
                            Sfr(Ev,H,A,S,P)=nanstd(fr2{Ev,A}*1000);
                            end
                            
                        else
                            Mfr{neuro}(Ev,H,A,S,P)=NaN;
                            Sfr(Ev,H,A,S,P)=NaN;
                        end
                    end
                    
                      if H==2 && S==1
                                clear temperar
                                j=1;
                                for i=1:4:length(Evidences)
                                    temperar{1,j}=fr{i,A};
                                    temperar{2,j}=Evidences(i);
                                    j=j+1;
                                end
                                Comparison{neuro,P,A}=temperar;
                                
                            end
                    
                    if~ Alice && S==1 %%ONLY PLOTTING CORRECT
                        nexttile(P);
                        plot(Evidences,Mfr{neuro}(:,H,A,S,P),col{H,S+2*(A-1)},'Markersize',5,'MarkerFaceColor',col2{H,S+2*(A-1)});
                        ylabel('FR (spikes/sec)')
                        xlabel('Sample Angle (deg)')
                        hold on
                        title(titly,'Interpreter', 'none');
                    else
                       
                      
                    end
                    
                end
                
            end
          disp('test')  
        end
        if ~isempty(fr{5,1}) && ~isempty(fr{5,2})
        pval(neuro,H)=ranksum(fr{5,1},fr{5,2});
        else
            pval(neuro,H)=nan;
        end
        
    end
    
    %Look for effect of Hrate*Active, acounting for effects of Ev, H, and
    %Active individually.
    for P=1:5
        T=[0,0,0,0;1,0,0,0;0,1,0,0; 0,0,1,0;0,1,1,0];fm=fitlm(X(X(:,4)==P,1:3),y(X(:,4)==P),T,'VarNames',{'Ev','Haz','Active','FR'});
        HrateEffect{neuro,P}=fm.Coefficients;
     
       
        for H=1:2
            for A=1:2
                miniData=[X(X(:,4)==P & X(:,2)==H & X(:,3)==A,1),y(X(:,4)==P & X(:,2)==H & X(:,3)==A,1)];
                HolderX=[];
            for Booties=1:100
                for Ev=1:length(Evidences)
                    if sum(miniData(:,1)==Evidences(Ev))>1
                    HolderX=[HolderX;datasample(miniData(miniData(:,1)==Evidences(Ev),:),sum(miniData(:,1)==Evidences(Ev)))];
                    else
                         HolderX=[HolderX;miniData(miniData(:,1)==Evidences(Ev),:)];
                    end
                end
               temp=polyfit(HolderX(:,1),HolderX(:,2),1);
               intercepts(Booties,A,H)=temp(2);
               slopes(Booties,A,H)=temp(1);
            end
            
                tempe=polyfit(miniData(:,1),miniData(:,2),1);
                RealInt(A,H)=tempe(2);
                RealSlope(A,H)=tempe(1);
                
                 LinerThings={'-r','-b','-m','-c'};
        nexttile(P);
        hold on
        plot(Evidences,RealInt(A,H)+Evidences*RealSlope(A,H),LinerThings{A+2*(H-1)});
      
           
            end
        end
   

        
     SignRankLowInt(neuro,P)=signrank((intercepts(:,1,1)-intercepts(:,2,1)));
     SignRankLowSlope(neuro,P)=signrank((slopes(:,1,1)-slopes(:,2,1)));
     SignRankNInt(neuro,P)=signrank((intercepts(:,1,2)-intercepts(:,2,2)));
     SignRankNSlope(neuro,P)=signrank((slopes(:,1,2)-slopes(:,2,2)));
     IntRealLowDiff(neuro,P)=RealInt(1,1)-RealInt(2,1);
     IntRealNDiff(neuro,P)=RealInt(1,2)-RealInt(2,2);
     SlpeRealLowDiff(neuro,P)=RealSlope(1,1)-RealSlope(2,1);
     SlpeRealNDiff(neuro,P)=RealSlope(1,2)-RealSlope(2,2);
     
       
    end
    
    
    if ~Alice
        %legend(['H=.05, A1; ', namer{1}], ['H=.05, A1; ', namer{2}], ['H=.05, A2; ', namer{1}], ['H=.05, A2; ', namer{2}],['H=.5, A1; ', namer{1}], ['H=.5, A1; ', namer{2}], ['H=.5, A2; ', namer{1}], ['H=.5, A2; ', namer{2}], 'Location', 'eastoutside')
           legend(['H=.05, A1; ', namer{1}], ['H=.05, A2; ', namer{1}],['H=.5, A1; ', namer{1}],  ['H=.5, A2; ', namer{1}], 'Location', 'eastoutside')

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

