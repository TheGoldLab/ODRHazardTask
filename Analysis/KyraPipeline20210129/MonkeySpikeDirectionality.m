function [ Fits, NumTri ] = MonkeySpikeDirectionality(data, figgy,b)
% reparse=1;
% if reparse
% ResponseCutoff=30; %(responses within this many degree of target are considered a response, rest are NaN)
% datahold=KAS_Parser_794_DelayNeural(ResponseCutoff,1);% KAS_Parser_794_recoded(ResponseCutoff);
% end
% clear x
%%
datahold=data;
FRlist=[];
Dirlist=[];
wrap=1;
if wrap
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
end

spikelist={};
for i=1:unique(datahold.Num_Neuron)
    name=['Neuron_',num2str(i)];
    spikelist=[spikelist,name];
end

Evidences=unique(datahold.Sample_Angle);
Evidences=Evidences(~isnan(Evidences));
NumTri=[];
for neuro=1:unique(datahold.Num_Neuron);
spiketimes=datahold{:,23+neuro};
clear y
for period=1:3
for Ev=1:length(Evidences) 

fr{Ev}=[];
Mfr=[];
Sfr=[];

end

num_sessions = size(unique(datahold{:,1}),1);

fileNAMES=unique(datahold{:,1});
for Sess=1:num_sessions
if period==1
    starter=datahold.SampleOn;%datahold.SacOn; %datahold.SampleOff;
    ender=datahold.SampleOff; %starter+500;%datahold.FPOff;
    titly={[fileNAMES{Sess}],['Visual'],['Neuron ',[datahold.Properties.VariableNames{23+neuro}]]};
elseif period==2
    starter=datahold.SampleOff;%datahold.SacOn; %datahold.SampleOff;
    ender=datahold.FPOff; %starter+500;%datahold.FPOff;
        titly={[fileNAMES{Sess}],['Memory'],['Neuron ',[datahold.Properties.VariableNames{23+neuro}]]};
else
    starter=datahold.SacOn; %datahold.SampleOff;
    ender=starter+500;%datahold.FPOff;
        titly={[fileNAMES{Sess}],['Motor'],['Neuron ',[datahold.Properties.VariableNames{23+neuro}]];};
end
  FRlist=[];
Dirlist=[];  
Baselist=[];


    
for Ev=1:length(Evidences) 

fr{Ev}=[];
Mfr=[];
Sfr=[];

end
             sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});

             if sum(cell2mat(spiketimes(sessionBreakdown)))
               for Ev=1:length(Evidences)
                TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev));
 NumTri(Ev,2)=length(TOI);
 NumTri(Ev,1)=Evidences(Ev);
                     for trials=1:length(TOI)
                         spks=sum([spiketimes{TOI(trials)}]>starter(TOI(trials)) & [spiketimes{TOI(trials)}]<ender(TOI(trials)));
                          fr{Ev}=[fr{Ev};spks/(ender(TOI(trials))-starter(TOI(trials)))];
                            Baselist=[Baselist,sum([spiketimes{TOI(trials)}]<[datahold.SampleOn(TOI(trials))])/[datahold.SampleOn(TOI(trials))]];
                     end
                     FRlist=[FRlist,fr{Ev}'];
                     Dirlist=[Dirlist,ones(size(fr{Ev}'))*Evidences(Ev)];
                    
                 Mfr(Ev,Sess)=nanmean(fr{Ev}*1000);
                Sfr(Ev,Sess)=nanstd(fr{Ev}*1000);
               end
               Baselist=Baselist*1000;
                FRlist=FRlist*1000;
                         figr=figure(figgy+neuro-1);
                         figr.Position=[680 558 1006 420];
                         subplot(1,3,period)
                         % x3= stretch factor, x2=mu, x1=kappa
                          Baseline=nanmean(Baselist); %This should eventually be calculated, but need to ask about ITI spikes first, for now is Baseline
%                              yfitM=@(x) x(3)*(exp(x(1)*cos(Evidences*pi/180-x(2))))./(2*pi*besseli(0,x(1)))+Baseline;    
%                              VonFitM = @(x)nansum((Mfr(:,Sess)-yfitM(x)).^2);
%                             % options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% %                              x(Sess,:)=fmincon(VonFitM,[sigma2kappa(15),pi/180*45,10],[],[],[],[],[0,0,0],[100,2*(pi),30]);
%                           
%                               yfit=@(x) x(3)*(exp(x(1)*cos(Dirlist*pi/180-x(2))))./(2*pi*besseli(0,x(1)))+Baseline;    
%                            VonFit = @(x)nansum((FRlist-yfit(x)).^2);
%                               y(period,:)=fmincon(VonFit,[sigma2kappa(15),pi/180*45,10,10],[],[],[],[],[1,0,0,0],[30,2*(pi),50,50]);
                         polar(Evidences*pi/180,Mfr(:,Sess));  
                         hold on
                        %  polar(Evidences*pi/180,yfitM(x(Sess,:)));
                         title(titly,'Interpreter', 'none');
                
 test(period,:)=VMFitting(Evidences,Mfr(:,Sess),Baseline);
                         Extend=[0:15:360];

                         VMFit=@(x) x(3)*(exp(x(1)*cos(Extend*pi/180-x(2))))./(2*pi*besseli(0,x(1)))+Baseline;
%                           polar(Extend*pi/180,VMFit(y(period,:)))
                          polar(Extend*pi/180,VMFit(test(period,:)))
                          %polar(Extend*pi/180,VMFit(x(Sess,:)))
%             x(Sess,2)=mod(x(Sess,2)/pi*180, 360);
%             x(Sess,1)=kappa2sigma(x(Sess,1));
            
%         y(period,2)=mod(y(period,2)/pi*180, 360);
%             y(period,1)=kappa2sigma(y(period,1));
%             y(period,4)=Baseline;
            
                    test(period,2)=mod(test(period,2)/pi*180, 360);
            test(period,1)=kappa2sigma(test(period,1));
          
            
            
             end
end
end
y=test;
y=num2cell(y);
 
 y=cell2table(y,'VariableNames',{'STD','Mean','Height','Baseline'},'RowNames',{'Visual', 'Memory','Motor'});
Fits{neuro}=y;

 print([b,'_neuron_', num2str(neuro), '_PETH_Directionality'],'-dpdf')
end
if ~isempty(NumTri)
NumTri=num2cell(NumTri);
NumTri=cell2table(NumTri, 'VariableNames', {'Direction', 'NumTrials'});
 Fits=cell2table(Fits,'VariableNames',spikelist)
else
    Fits=[];
    
end

%%
% figure
% plot(Evidences*pi/180, Mfr(:,Sess));
% hold on
% plot(Extend*pi/180,VMFit([sigma2kappa(x(1)),pi/180*x(2),x(3)]));%[sigma2kappa(30),pi/180*45,20]