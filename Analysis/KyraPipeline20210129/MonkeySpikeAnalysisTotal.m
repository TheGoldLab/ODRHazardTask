function [  ] = MonkeySpikeAnalysisTotal(data, figgy, b)
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
col={'r','b'};

%First, need to split up by evidence, maybe later also by if chose A or B
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
Evidences=unique(ceil(datahold.Sample_Angle));
datahold.Sample_Angle=ceil(datahold.Sample_Angle);
Evidences=Evidences(~isnan(Evidences));
Choices=unique(datahold.Choice);
Choices=Choices(~isnan(Choices));
colors={'r','--r','b','--b'};

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


for Sess=1%:num_sessions
            %Find out which entries in the mega file are from a particular
            %session
         sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
         if sum(cell2mat(spiketimes(sessionBreakdown)))
for H=1:2;
    liner=[];
                 for A=1:2;
     
          for Ev=1:length(Evidences)
              
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
          
              TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev) & abs(datahold.H_Rate-Hs(H))<.0001 & subsettle==1);
        % figure(A)
%               subplot(5,2,2*Ev-(2-H))
%               title({['Evidence= ', num2str(Evidences(Ev))]; ['H=  ',num2str(Hs(H))]});
%               hold on
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
                 bins=[-500:100:3000];
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
          Wantpsth=1;
          if Wantpsth
                      fig=figure(figgy+neuro-1);
                      fig.Position=[680 558 560 1000];
                      hold on
               subplot(5,2,2*Ev-(2-H))
               MNSpkCnt=nanmean(PSTH{Ev,H},1)/200*1000;
               baseline=0; %nanmean(MNSpkCnt(1:4));
               hold on
               if ~isempty(MNSpkCnt)
               liner(A)=plot(bins(2:end-1),MNSpkCnt-baseline, col{A});
               else
                  liner(A)=plot(0,0, col{A}); 
               end
         
             title({['Evidence= ', num2str(Evidences(Ev))]; ['H=  ',num2str(Hs(H))];['Neuron ' ,[datahold.Properties.VariableNames{23+neuro}]]},'Interpreter', 'none');
if H==1 && 2*Ev-(2-H)==1
             legend('Active=45', 'Active=315','location','southwest');
elseif H==2 && 2*Ev-(2-H)==2
             legend('Prev Active=45', 'Prev Active=315','location','southwest');

end
if A==2
plot([1532,1532],[-1,20],'--k');
end
%           end
              end  
          end
          ylabel('Avg FR');
          xlabel('Time(ms) From Samp Off');
          
          axis([-1000,3000,-10,70]);
             
            for Ev=1:length(Evidences) 
                Mfr(Ev,H,A)=nanmean(fr{Ev,H}*1000);
                Sfr(Ev,H,A)=nanstd(fr{Ev,H}*1000);
            end
            wantSummary=0;
            if wantSummary
            figure(5) 
            hold on
             errorbar(Evidences,Mfr(:,H,A),Sfr(:,H,A),colors{H+2*(A-1)});
             legend('H=.05,Active=45', 'H=.05, Active=315', 'H=.95, Prev Active=45', 'H=.95 Prev Active=315');
             xlabel('Evidence');
             ylabel('FR(Sp/s)');
            end
             end
             end
             
         end
end    
set(gcf, 'PaperPositionMode', 'auto');
 print([b,'_neuron_', num2str(neuro), '_PETH_Full'],'-dpdf')
end
end
