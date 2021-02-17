function [ NumTri ] = MonkeySpikeAnalysisSubset(data, subby, figgy,b)
%Pick/parse actual data files
% reparse=1;
% if reparse
% ResponseCutoff=30; %(responses within this many degree of target are considered a response, rest are NaN)
% datahold=KAS_Parser_794_DelayNeural(ResponseCutoff,2);% KAS_Parser_794_recoded(ResponseCutoff);
% end
%%
% Goals:  Make raster plot of spike times, make PTSH, Look at tuning (avg
% normalized H rate)
datahold=data;
clear Mfr
clear Sfr
col={'r','--r','b','--b'};

%First, need to split up by evidence, maybe later also by if chose A or B
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
Evidences=unique(ceil(datahold.Sample_Angle));
datahold.Sample_Angle=ceil(datahold.Sample_Angle);
Evidences=Evidences(~isnan(Evidences));
Choices=unique(datahold.Choice);
Choices=Choices(~isnan(Choices));

NumTri=[];
clear fr
num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});
allFR=cell(1,2);
for neuro=1:unique(datahold.Num_Neuron);
spiketimes=datahold{:,23+neuro};
for Ev=1:length(Evidences) 
    for H=1:2
fr{Ev,H}=[];
Rastas{Ev}=[];
PSTH{Ev,H}=[];
    end
end

Hs=[.05,.95];
wantChoices=1;


if subby==1 %corr vs err
    opts=[0,1];
    catagory=datahold.Correct;
    namer={'Err','Corr'};
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

for Sess=1%:num_sessions
            %Find out which entries in the mega file are from a particular
            %session
         sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
         if sum(cell2mat(spiketimes(sessionBreakdown)))
for H=1:2;
   
     
     liner=[];
          for Ev=1:length(Evidences)
              for A=1:2;
                for S=1:2; 
              if wantChoices
                  if H==1
                  subsettle= datahold.Active_Angle==Choices(A) & datahold.TACP_Actual>5; % Active_Angle
                  else
                      prevAngle=datahold.Active_Angle;
                      prevAngle=[NaN;prevAngle(1:end-1)];
                   subsettle= prevAngle==Choices(A); % Active_Angle

                  end
              else
                  subsettle=ones(size(datahold.Sample_Angle));
                  A=1;
              end
          
              TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1 & catagory==opts(S));
       NumTri(Ev,(4*(H-1)+2*(A-1)+S))=length(TOI);
              
      
              % figure(A)
%               subplot(5,2,2*Ev-(2-H))
%               title({['Evidence= ', num2str(Evidences(Ev))]; ['H=  ',num2str(Hs(H))]});
%               hold on
   bins=[-500:100:3000];
              for trials=1:length(TOI)

              %Raster
            
                 spiky=spiketimes{TOI(trials)};
                 spiky=spiky(spiky<4000);
                 wantRasta=0;
              if datahold.SampleOn(TOI(trials))<5000 & datahold.SampleOff(TOI(trials))<5000 & datahold.FPOff(TOI(trials)) & wantRasta
                  scatter(spiky,trials*ones(1,length(spiky)),'.k');
                  scatter(datahold.SampleOn(TOI(trials)),trials,'db');
                  scatter(datahold.SampleOff(TOI(trials)),trials,'dg');
                  scatter(datahold.FPOff(TOI(trials)),trials,'dr');
              end

              
              %PSTH
                 spiky= spiketimes{TOI(trials)};
                 spiky=spiky(spiky<4000)-datahold.SampleOff(TOI(trials));
              
                 binnedspikes=[];
                 for bin=1:length(bins)-2
                     binnedspikes(1,bin)=sum(spiky>bins(bin) & spiky<bins(bin+2));
                 end
                 
                       PSTH{Ev,H}=[PSTH{Ev,H};binnedspikes];
                            
           
       

              
              %Avg FR 
              spks=sum([spiketimes{TOI(trials)}]>datahold.SampleOff(TOI(trials)) & [spiketimes{TOI(trials)}]<datahold.FPOff(TOI(trials)));
              fr{Ev,H}=[fr{Ev,H};spks/(datahold.FPOff(TOI(trials))-datahold.SampleOff(TOI(trials)))];
              allFR{A}=[allFR{A};spks/(datahold.FPOff(TOI(trials))-datahold.SampleOff(TOI(trials)))];
              end
%               if ~isempty(TOI)
          ylabel('Trial');
          xlabel('time(ms)');
          axis([-1000,3000,-10,70]);
          Wantpsth=1;
          if Wantpsth
                      fig=figure(figgy+neuro-1);
                      fig.Position=[680 558 560 1000];
                      hold on
               subplot(5,2,2*Ev-(2-H))
               if ~isempty(PSTH{Ev,H})
               MNSpkCnt=nanmean(PSTH{Ev,H},1)/200*1000;
               else
                   MNSpkCnt=nans(1,length(bins)-2);
               end
               
               PSTH{Ev,H}=[];
               baseline=0; %nanmean(MNSpkCnt(1:4));
               hold on
             liner(A+2*(S-1))= plot(bins(2:end-1),MNSpkCnt-baseline, col{S+2*(A-1)});
         
             title({['Evidence= ', num2str(Evidences(Ev))]; ['H=  ',num2str(Hs(H))];['Neuron ' datahold.Properties.VariableNames{23+neuro}]},'Interpreter', 'none');
if H==1 && 2*Ev-(2-H)==1 && A+2*(S-1)==4
            legend(['Active=45 ',namer{1}], ['Active=45 ',namer{2}],['Active=315 ',namer{1}], ['Active=315 ',namer{2}],'location','southwest');
elseif H==2 && 2*Ev-(2-H)==2 
             legend(['Prev Active=45 ',namer{1}], ['Prev Active=45 ',namer{2}],['Prev Active=315 ',namer{1}], ['Prev Active=315 ',namer{2}],'location','southwest');

end
if A==2
plot([1532,1532],[-1,20],'--k');
end
%           end
              end  
          end
          ylabel('Avg FR');
          xlabel('Time(ms) From Samp Off');
%              
%             for Evs=1:length(Evidences) 
%                 Mfr(Evs,H,A)=nanmean(fr{Evs,H}*1000);
%                 Sfr(Evs,H,A)=nanstd(fr{Evs,H}*1000);
%             end
%             wantSummary=0;
%             if wantSummary
%             figure(5) 
%             hold on
%              errorbar(Evidences,Mfr(:,H,A),Sfr(:,H,A),colors{H+2*(A-1)});
%              legend('H=.05,Active=45', 'H=.05, Active=315', 'H=.95, Prev Active=45', 'H=.95 Prev Active=315');
%              xlabel('Evidence');
%              ylabel('FR(Sp/s)');
%             end
        end
    end
             end
             
         end
        
end
set(gcf, 'PaperPositionMode', 'auto');

 print([b,'_neuron_', num2str(neuro), '_PETH_',label],'-dpdf')
if neuro==datahold.Num_Neuron
NumTri=num2cell(NumTri);
NumTri=cell2table(NumTri,'VariableNames',{['HL_A1_',namer{1}],['HL_A1_',namer{2}],['HL_A2_',namer{1}],['HL_A2_',namer{2}],['HH_A1_',namer{1}],['HH_A1_',namer{2}],['HH_A2_',namer{1}],['HH_A2_',namer{2}]},'Rownames',{['Ev_',num2str(Evidences(1))],['Ev_',num2str(Evidences(2))],['Ev_',num2str(Evidences(3))],['Ev_',num2str(Evidences(4))],['Ev_',num2str(Evidences(5))]})
end
end
end
