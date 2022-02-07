function [ Fits, MFRTotal, SFRTotal, NumTri, pval,Index, Comparison ] = MonkeySpikeDirectionality(data, figgy,b)
% reparse=1;
% if reparse
% ResponseCutoff=30; %(responses within this many degree of target are considered a response, rest are NaN)
% datahold=KAS_Parser_794_DelayNeural(ResponseCutoff,1);% KAS_Parser_794_recoded(ResponseCutoff);
% end
% clear x
%%
datahold=data;
FRlist=[];
FRbaseList=[];
Dirlist=[];


%Name of all the neurons for ease of titles later
spikelist={};
for i=1:unique(datahold.Num_Neuron)
    name=['Neuron_',num2str(i)];
    spikelist=[spikelist,name];
end

Evidences=unique(datahold.Sample_Angle);
Evidences=Evidences(~isnan(Evidences));
NumTri=[];
Sess=1;
for neuro=1:unique(datahold.Num_Neuron);

    spiketimes=datahold{:,24+neuro};
    clear y
    figr=figure(figgy+neuro-1);
    tiledlayout(2,4, 'TileSpacing', 'compact')    ;
    for P=1:4
        FRlist=[];
        FRbaseList=[];
        Dirlist=[];
        num_sessions = size(unique(datahold{:,1}),1);
        fileNAMES=unique(datahold{:,1});
        if P==1
            starter=datahold.SampleOn;%datahold.SacOn; %datahold.SampleOff;
            ender=starter+50;%datahold.FPOff;
            titly={[fileNAMES{Sess}],['Visual'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]]};
        elseif P==2
            starter=datahold.SampleOn+50; %datahold.SampleOff;%datahold.SacOn; %datahold.SampleOff;
            ender=datahold.FPOff; %starter+500;%datahold.FPOff;
            titly={[fileNAMES{Sess}],['Memory'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]]};
        elseif P==3
            starter=datahold.FPOff; %datahold.SampleOff;
            ender=datahold.FPOff+500;datahold.TargAcq;%datahold.FPOff;
            titly={[fileNAMES{Sess}],['Motor'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]];};
            
        else
            starter=datahold.FPOff+500;%datahold.TargAcq; %datahold.SampleOff;
            ender=datahold.FPOff+1000; %datahold.TargAcq+500;%datahold.FPOff;
            titly={[fileNAMES{Sess}],['Reward'],['Neuron ',[datahold.Properties.VariableNames{24+neuro}]];};
        end
        
        basestart=datahold.SampleOn-500;
        basestop=datahold.SampleOn;
        
        for Ev=1:length(Evidences)
            
            fr{Ev}=[];
            frbase{Ev}=[];
            
            sessionBreakdown=strcmp(datahold{:,1},fileNAMES{Sess});
            
            TOI=find(sessionBreakdown & datahold.Sample_Angle==Evidences(Ev));
            NumTri(Ev,2)=length(TOI);
            NumTri(Ev,1)=Evidences(Ev);
            for trials=1:length(TOI)
                spks=sum([spiketimes{TOI(trials)}]>starter(TOI(trials)) & [spiketimes{TOI(trials)}]<ender(TOI(trials)));
                fr{Ev}=[fr{Ev};spks/(ender(TOI(trials))-starter(TOI(trials)))];
                spks=sum([spiketimes{TOI(trials)}]>basestart(TOI(trials)) & [spiketimes{TOI(trials)}]<basestop(TOI(trials)));
                frbase{Ev}=[frbase{Ev};spks/(basestop(TOI(trials))-basestart(TOI(trials)))];

            end
            FRlist=[FRlist,fr{Ev}'];
            FRbaseList=[FRbaseList,frbase{Ev}'];
            Dirlist=[Dirlist,ones(size(fr{Ev}'))*Evidences(Ev)];
            
            Mfr(Ev)=nanmean(fr{Ev}*1000);
            Sfr(Ev)=nanstd(fr{Ev}*1000)/sqrt(length(fr{Ev}));
        end
        clear temperar
        for i=1:length(Evidences)
        temperar{1,i}=fr{i};
        temperar{2,i}=Evidences(i);
        end
        Comparison{neuro,P}=temperar;
        FRlist=FRlist*1000;
                FRbaseList=FRbaseList*1000;

        Baseline=0;
        
        [pval(neuro,P),~]=anova1(FRlist,Dirlist,'off')
        
        %Determine if directional
        
        
        
        nexttile(P)
        polar(Evidences*pi/180,Mfr(:));
        hold on
        titly{4}=['pvalue= ', num2str(pval(neuro,P))];
        title(titly,'Interpreter', 'none');
        
        test(P,:)=VMFitting(Evidences,Mfr(:),Baseline);
        Extend=[0:15:360];
        
        VMFit=@(x) x(3)*(exp(x(1)*cos(Extend*pi/180-x(2))))./(2*pi*besseli(0,x(1)))+Baseline;
        polar(Extend*pi/180,VMFit(test(P,:)))
        
        %Non-polar versions for funzies
        nexttile(4+P)
        errorbar(Evidences,Mfr(:),Sfr(:));
        hold on
        title(titly,'Interpreter', 'none');
        Index(neuro,P)=(max(VMFit(test(P,:)))-min(VMFit(test(P,:))))/(max(VMFit(test(P,:)))+min(VMFit(test(P,:))))
        
        plot(Extend,VMFit(test(P,:)));
        ylabel('FR')
        xlabel('Target Angle (deg)')
        xlim([0,360]);
        
        %Convert to degrees and std for ease of interpretation
        test(P,2)=mod(test(P,2)/pi*180, 360);
        test(P,1)=kappa2sigma(test(P,1));
        MFRTotal(P,:,neuro)=Mfr;
        SFRTotal(P,:,neuro)=Sfr;
    end
    y=test;
    y=num2cell(y);
    
    y=cell2table(y,'VariableNames',{'STD','Mean','Height','Baseline'},'RowNames',{'Visual', 'Memory','Motor','Reward'});
    Fits{neuro}=y;
    
 set(gcf, 'Units', 'Inches')
     set(gcf,'Position',[ 2.4167    1.1583   10.9250    6.6500]);
pos=get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto','PaperUnits', 'Inches', 'PaperSize', [pos(3),pos(4)]);
saveas(gcf,[b,'_neuron_', num2str(neuro), '_Directionality.pdf'])  
end
if ~isempty(NumTri)
    NumTri=num2cell(NumTri);
    NumTri=cell2table(NumTri, 'VariableNames', {'Direction', 'NumTrials'});
    Fits=cell2table(Fits,'VariableNames',spikelist)
else
    Fits=[];
    
end
end
