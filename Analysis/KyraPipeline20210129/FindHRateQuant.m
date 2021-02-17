function [ z, NumTrialsLLR ] = FindHRateQuant(data,figgy)

%Pick/parse actual data files
% reparse=1;
% if reparse
% ResponseCutoff=30; %(responses within this many degree of target are considered a response, rest are NaN)
% datahold=KAS_Parser_794_DelayNeural(ResponseCutoff,2);% KAS_Parser_794_recoded(ResponseCutoff);
% end

clear NumTrialsLLR
datahold=data;
%% General conditions
wantChoiceOn=0; %0= all, 1= correct answer on, 2=correct answer off during sample period
STRENGTH=0;
wantExcludeBegining=0;
basedOnALLTrials=0;
Cutoff=150;
datahold{:,9}=ceil(datahold{:,9});
datahold{:,16}=ceil(datahold{:,16});

LLRS=unique(datahold{:,9});
adjustforzeroH=(unique(datahold{:,16}));
adjustforzeroH=adjustforzeroH(~isnan(unique(datahold{:,16})));
QuantConv= [log(1/54),log(10/20),log(15/15),log(20/10),log(54/1)]; %[-3.277    -1.3863         0   1.3863   3.277];%[-2.3026   -0.6931         0    0.6931    2.3026];
for i=1:length(LLRS)
datahold{datahold{:,9}==LLRS(i),9}=QuantConv(i);
datahold{datahold{:,16}==adjustforzeroH(i),16}=QuantConv(i);
end
%%
if wantChoiceOn==1
colors={'ro','m-','bo','c-','b','k'};

elseif wantChoiceOn==2
 colors={'g','c','b','k','b','k'};  
else
    colors={'ro','m-','bo','c-','b','k'};
end

num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});

%General Settings:
Hs=unique(datahold{:,3});%Hs=[.05,.5,.5,.9];
ChoiceNoiseFlag=1;
% CalcH=nans(4,length(breakdown));

clear z; clear x; clear y;j=1;

%Session Breakdown condition 
wantSessions=1;
if wantSessions
breakdown=[1:num_sessions];
else
    breakdown=num_sessions+1;
end    



%Learning exclusion 
if wantExcludeBegining
    if basedOnALLTrials
            %Based on trials within session
                 learningexclusion=datahold.Trial_Within_Session>Cutoff;
    else
         %if want based on by unbroken trials
         [~,B,~]=unique(datahold.Session);
         learningexclusion=ones(size(datahold.Session));
    
            for k=1:length(B)
            learningexclusion(B(k):B(k)+Cutoff)=0;
            end
    
        %For if change H rate in session
        holdingtemp=diff(datahold.H_Rate);
        switchedAng=find(holdingtemp>0);
        for k=1:length(switchedAng)
        leanringexclusion(switchedAng(k):switchedAng(k)+Cutoff)=0;
        end
    end
else
    learningexclusion=ones(size(datahold.Trial_Within_Session));
end





%% 


   for Session=1:length(breakdown)
    firstj=j;
for Hrate=1:length(unique(datahold.H_Rate))
%
pChange=[];

  if breakdown==num_sessions+1
    sessionBreakdown=ones(size(datahold.Noise));
    name='All';
  else
     sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Session});
    name=fileNAMES{Session};
  end

%Parse when picked Targ 1 vs Targ 2 for fitting
ChooseTarg1=[];
ofInt=sessionBreakdown & learningexclusion  &  abs(datahold.H_Rate-Hs(Hrate))<.0000000000001 & abs(datahold.LLR)<4.5;

%Work only with session(s) of interest

datahold2=datahold(ofInt,:);
ChooseTarg1(datahold2.Choice==datahold2.T1_Angle)=1;
ChooseTarg1(isnan(datahold2.Choice))=nan;
ChooseTarg1(datahold2.Choice==datahold2.T2_Angle)=0;
ChooseTarg1=ChooseTarg1';
Actualchoices=~isnan(ChooseTarg1);
LLRs=datahold2.LLR(Actualchoices);
FullLLRs=datahold2.LLR;
GroundTruth=datahold2.Active_Angle==datahold2.T1_Angle;
Previous=[0;GroundTruth(1:end-1)];

    LLRForChange=datahold2.LLR_for_Change;
    Switch=[];
    Switch(:,j)=datahold2.Switch;
    ChooseT1=[];
    if ~isempty(ChooseTarg1)
    ChooseT1(:,j)=ChooseTarg1;
    end
    
    
%Figure out how many trials you have to work with:
datahold2.Sample_Angle(datahold2.Sample_Angle>180)=datahold2.Sample_Angle(datahold2.Sample_Angle>180)-360;
Evidences=unique(round(datahold2.Sample_Angle,-1));


Evidences=Evidences(~isnan(Evidences));
Choices=unique(datahold2.Choice);
Choices=Choices(~isnan(Choices));
for Active=1:2
for i=1:length(Evidences)
if Hrate==1
    NumTrialsLLR(Active+2*(Hrate-1),i, Session)=sum(datahold2.Active_Angle==Choices(Active) & round(datahold2.Sample_Angle,-1)==Evidences(i));
else
     prevAngle=datahold2.Active_Angle;
                      prevAngle=[NaN;prevAngle(1:end-1)];

   NumTrialsLLR(Active+2*(Hrate-1),i, Session)=sum(prevAngle==Choices(Active) & round(datahold2.Sample_Angle,-1)==Evidences(i));

end

end
end

%% Start to make scatter plot

if ~isempty(Actualchoices)

%     strong_evidence_before=find(abs(FullLLRs)>STRENGTH);
PreviousCorrect=0;
if PreviousCorrect
strong_evidence_before=find(abs(datahold2.LLR)>=STRENGTH & datahold2.Correct==1);
else 
strong_evidence_before=find(abs(datahold2.LLR)>=STRENGTH);
end
    strong_evidence_before=strong_evidence_before+1;
 
    
 if wantChoiceOn==1
     Strong_and_vis=intersect(find(datahold2.Correct_On(:)==1),strong_evidence_before);
elseif wantChoiceOn==2
Strong_and_vis=intersect(find(datahold2.Correct_On(:)==0),strong_evidence_before);
 else
     Strong_and_vis=strong_evidence_before;
 end    

           strong_evidence_before=intersect(Strong_and_vis,find(Actualchoices));

    disp('length')
    length(strong_evidence_before)

    %Scatter Plot
    figure(figgy)%(2*(Session-1)+1)
    xlabel('LLR for Change');
    ylabel('Probability of Switch');
    hold on;

%     scatter(LLRForChange(strong_evidence_before),Switch(strong_evidence_before,j));
%     title({['Actual H=', num2str(H),  ' NoiseChoice = ' num2str(NoiseProb), ' NoiseP= ' num2str(NoiseProb)]})

%% Begin averaging in bins to make psychometric function
pChange=[];

bins=unique(datahold2.LLR);
    timesResample=1000;

    
    for i=1:length(bins)
        pChange(i)=nanmean(Switch(strong_evidence_before(LLRForChange(strong_evidence_before)==bins(i)))); 
        errbarholder(i)=bootstrapCurve( timesResample,Switch(strong_evidence_before(LLRForChange(strong_evidence_before)==bins(i))));

    end
  %%  
    if nnz(isnan(pChange(:)))>0 
    nangone=find(~isnan(pChange(:)));
    pChangetemp=pChange(nangone);
    errbartemp=errbarholder(nangone);
    bins=bins(nangone);
    else
        pChangetemp=pChange;
       errbartemp=errbarholder;

    end
    psychofit=1;
    if psychofit
  %Psychometric curve
    figure(figgy) %(2*(Session-1)+1)
%     title({['Model threshold=   ', num2str(x(Session,2))]})%['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)], 

    errorbar(bins,pChangetemp, errbartemp, colors{1+2*(Hrate-1)}, 'Linewidth', 3);
    hold on

   
    yfit=@(x) x(1)+(1-2*x(1))./(1+exp((bins'-x(2)).*-x(3)));    
    LogisticFit = @(x)nansum((pChangetemp-yfit(x)).^2);
    x(Session,:,Hrate)=fmincon(LogisticFit,[0,.5,0],[],[],[],[],[0,-4,0],[.3,4,10]);
    %Showfit
    tempper=x(Session,:,Hrate);
    
    showfitX=-4:.2:4;
    showfitY=tempper(1)+(1-2*tempper(1))./(1+exp((showfitX-tempper(2)).*-tempper(3)));
    
    plot(showfitX,showfitY,  colors{2+2*(Hrate-1)}, 'Linewidth', 2');
        axis([-4,4,0,1]);
       
% 
%       plot(bins,yfit(x(Session,:,Hrate)), colors{2+2*(Hrate-1)});
if Hrate==2
       legend('Data, H=.05','Logistic Fit, H=.05' , 'Data, H=.5', 'Logistic Fit, H=.5','location','northwest');
                  set(gca,'FontSize',20)
else
           legend(['Data, H=', num2str(Hs(Hrate))],['Logistic Fit, H=' ,num2str(Hs(Hrate))], 'Data, H=.5', 'Logistic Fit, H=.5','location','northwest');
                  set(gca,'FontSize',20)
end

    end

%% Fit Model based on the H rate that maximizes the probability that the actual choices were chosen, semi-conditionalized on whether the correct target was shown or not


      if length(Actualchoices)>50
 %%%***** This version gives only a continuous list of the trials in which answer was shown or not         
% if wantChoiceOn==1
% Actualchoices=Actualchoices & datahold2.Correct_On==1;
% elseif wantChoiceOn==2
% Actualchoices=Actualchoices & datahold2.Correct_On==0;
% 
% end   

%This version puts nans as choices if the answer was given so no weight in
%the minimization 
if wantChoiceOn==1
ChooseT1(datahold2.Correct_On==0,j)=nan;
elseif wantChoiceOn==2
ChooseT1(datahold2.Correct_On==1,j)=nan;

end   


sum(Actualchoices)
z(:,Hrate)=MSFitting(FullLLRs(Actualchoices),ChooseT1(Actualchoices,j));
        else
          z(:,Hrate)=[nan;nan;nan];
end
      end

      
      
% %% Plot error landscape for model fits      
%  prob3=@(x) PredictedProb(FullLLRs(Actualchoices),x(1),x(2),x(3));
% 
% error3= @(x) -nansum((1-ChooseT1(Actualchoices,j)).*log(1.0000-prob3(x))+ChooseT1(Actualchoices,j).*log(prob3(x)+.00000)); %Equation 15
% 
% 
% 
% figure(3)
% teting=[]; 
% for i=.05:.01:.5
% teting=[teting,error3([i,z(2,Session,Hrate), z(3,Session,Hrate)])];
% end
% plot([.05:.01:.5],teting)%,colors{j}
%  hold on
% teting=[];
% for i=.05:.05:1
% teting=[teting,error([i,x(2)*10,Noise*10])];
% end
% plot([.05:.05:1],teting, 'r')
% 
% xlabel('H');
% % title({ ['Model H=', num2str(z(1,Session,Hrate))]}) %['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)],
% j=j+1;
    
end
             figure(figgy)%(2*(Session-1)+1)
             if Hrate==2
    title({['H=.05, Model Fit HRate=   ', num2str(z(1,1))],['H=.5, Model Fit HRate=   ', num2str(z(1,2))],[fileNAMES{Session}]},'Interpreter', 'none')%['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)], 
             else
                     title({['H=', num2str(Hs(Hrate)), ' Model Fit HRate=   ', num2str(z(1,1))],[fileNAMES{Session}]},'Interpreter', 'none')%['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)], 

             end
   end
 avgfits=nanmean(z,2);
 z;
 z=num2cell(z');
 NumTrialsLLR;
 NumTrialsLLR=num2cell(NumTrialsLLR);
 if length(Evidences)==5
 if size(NumTrialsLLR,1)==2 && Hs<=.5
      NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','M2','C','M1','T1'},'RowNames',{'T1_HL', 'T2_HL'});
 z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'Low_H'});
 elseif size(NumTrialsLLR,1)==2 && Hs>.5
      NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','M2','C','M1','T1'},'RowNames',{'T1_HH', 'T2_HH'});
    z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'High_H'});
 else
 NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','M2','C','M1','T1'},'RowNames',{'T1_HL', 'T2_HL', 'T1_HH', 'T2_HH'});
     z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'Low_H','High_H'});

 end
 end
 end