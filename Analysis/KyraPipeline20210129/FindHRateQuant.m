function [ z, NumTrialsLLR ] = FindHRateQuant(data,figgy,b)

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
datahold.Sample_Angle(datahold.Sample_Angle>180)=datahold.Sample_Angle(datahold.Sample_Angle>180)-360;
datahold{:,8}=round(datahold{:,8});
datahold{:,16}=abs(datahold{:,8}).*sign(datahold{:,16});

Locations=unique(datahold{:,8});
LocationsForChange=unique(datahold{:,16});
LocationsForChange=LocationsForChange(~isnan(LocationsForChange));
EvidenceForChange=datahold{:,16};
Evidence=datahold{:,8};
LLConv= [log(1/54),log(2/54),log(3/54),log(10/30),log(20/20),log(30/10),log(54/3),log(54/2),log(54/1)]; %[-3.277    -1.3863         0   1.3863   3.277];%[-2.3026   -0.6931         0    0.6931    2.3026];
for i=1:length(LocationsForChange)
    EvidenceForChange(EvidenceForChange==LocationsForChange(i))=LLConv(i);
    Evidence(Evidence==Locations(i))=LLConv(i);
end
%%

colors={'ro','m-','bo','c-','b','k'};

num_sessions = size(unique(datahold{:,1}),1);
fileNAMES=unique(datahold{:,1});

%General Settings:
Hs=unique(datahold{:,3});%Hs=[.05,.5,.5,.9];

clear z; clear x; clear y;j=1;


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
Session=1;
firstj=j;
for Hrate=1:length(unique(datahold.H_Rate))
    %
    pChange=[];
    
    
    sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Session});
    name=fileNAMES{Session};
    
    
    %Parse when picked Targ 1 vs Targ 2 for fitting
    ChooseTarg1=[];
    ofInt=sessionBreakdown & learningexclusion  &  abs(datahold.H_Rate-Hs(Hrate))<.00001 & abs(datahold.LLR)<4.5;

    
    datahold2=datahold(ofInt,:);
    ChooseTarg1(datahold2.Choice==datahold2.T1_Angle)=1;
    ChooseTarg1(isnan(datahold2.Choice))=nan;
    ChooseTarg1(datahold2.Choice==datahold2.T2_Angle)=0;
    ChooseTarg1=ChooseTarg1';
    Actualchoices=~isnan(ChooseTarg1);
    EFT2=EvidenceForChange(ofInt);
    LLRs=EFT2(Actualchoices );
    FullLLRs=Evidence(ofInt);
    GroundTruth=datahold2.Active_Angle==datahold2.T1_Angle;
    Previous=[0;GroundTruth(1:end-1)];
    
    Switch=[];
    Switch(:,j)=datahold2.Switch;
    ChooseT1=[];
    if ~isempty(ChooseTarg1)
        ChooseT1(:,j)=ChooseTarg1;
    end
    
    
    %Figure out how many trials you have to work with:
    datahold2.Sample_Angle(datahold2.Sample_Angle>180)=datahold2.Sample_Angle(datahold2.Sample_Angle>180)-360;
    Locations=unique(round(datahold2.Sample_Angle,-1));
    
    
    Locations=Locations(~isnan(Locations));
    Choices=unique(datahold2.Choice);
    Choices=Choices(~isnan(Choices));
    for Active=1:2
        for i=1:length(Locations)
            if Hrate==1
                NumTrialsLLR(Active+2*(Hrate-1),i, Session)=sum(datahold2.Active_Angle==Choices(Active) & round(datahold2.Sample_Angle,-1)==Locations(i));
            else
                ActiveAngle=datahold2.Active_Angle;
                NumTrialsLLR(Active+2*(Hrate-1),i, Session)=sum(ActiveAngle==Choices(Active) & round(datahold2.Sample_Angle,-1)==Locations(i));
                
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
        
        
        
        %% Begin averaging in bins to make psychometric function
        pChange=[];
        
        bins=LLConv;
        timesResample=1000;
        
        
        for i=1:length(bins)
            pChange(i)=nanmean(Switch(strong_evidence_before(EFT2(strong_evidence_before)==bins(i))));
            errbarholder(i)=bootstrapCurve( timesResample,Switch(strong_evidence_before(EFT2(strong_evidence_before)==bins(i))));
            
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
            LogisticFit = @(x)nansum((pChangetemp-yfit(x)').^2);
            x(Session,:,Hrate)=fmincon(LogisticFit,[0,.5,0],[],[],[],[],[0,-4,0],[.3,4,10]);
            %Showfit
            tempper=x(Session,:,Hrate);
            
            showfitX=-4:.2:4;
            showfitY=tempper(1)+(1-2*tempper(1))./(1+exp((showfitX-tempper(2)).*-tempper(3)));
            
            plot(showfitX,showfitY,  colors{2+2*(Hrate-1)}, 'Linewidth', 2');
            xlabel('LLR For Change');
            ylabel('P(Switch)');
         
            axis([-4,4,0,1]);
            
            %
            %       plot(bins,yfit(x(Session,:,Hrate)), colors{2+2*(Hrate-1)});
            if Hrate==2
                legend('Data, H=.05','Logistic Fit, H=.05' , 'Data, H=.5', 'Logistic Fit, H=.5','location','northwest');
                set(gca,'FontSize',14)
            else
                legend(['Data, H=', num2str(Hs(Hrate))],['Logistic Fit, H=' ,num2str(Hs(Hrate))], 'Data, H=.5', 'Logistic Fit, H=.5','location','northwest');
                set(gca,'FontSize',14)
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
    
    
end
figure(figgy)%(2*(Session-1)+1)
if Hrate==2
    title({['H=.05, Model Fit HRate=   ', num2str(z(1,1))],['H=.5, Model Fit HRate=   ', num2str(z(1,2))],[fileNAMES{Session}]},'Interpreter', 'none')%['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)],
else
    title({['H=', num2str(Hs(Hrate)), ' Model Fit HRate=   ', num2str(z(1,1))],[fileNAMES{Session}]},'Interpreter', 'none')%['Actual H=', num2str(H), ' NoiseChoice = ' num2str(NoiseProb), ' lapse= ' num2str(lapse)],
    
end

 set(gcf, 'Units', 'Inches')
 set(gcf,'Position',[7.8000    1.9167    6.8583    6.2250]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,[b,'_PsychometricFxn.pdf'])  


avgfits=nanmean(z,2);
z;
z=num2cell(z');
NumTrialsLLR;
NumTrialsLLR=num2cell(NumTrialsLLR);
if length(Locations)==9
    if size(NumTrialsLLR,1)==2 && Hs<=.5
        NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','T2_M','M2','T2_C','C','T1_C','M1', 'T1_M', 'T1'},'RowNames',{'T1_HL', 'T2_HL'});
        z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'Low_H'});
    elseif size(NumTrialsLLR,1)==2 && Hs>.5
        NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','T2_M','M2','T2_C','C','T1_C','M1', 'T1_M', 'T1'},'RowNames',{'T1_HH', 'T2_HH'});
        z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'High_H'});
    else
        NumTrialsLLR=cell2table(NumTrialsLLR,'VariableNames',{'T2','T2_M','M2','T2_C','C','T1_C','M1', 'T1_M', 'T1'},'RowNames',{'T1_HL', 'T2_HL', 'T1_HH', 'T2_HH'});
        z=cell2table(z,'VariableNames',{'Fit_H','Fit_Noise','Fit_Lapse'},'RowNames',{'Low_H','High_H'});
        
    end
end
end
