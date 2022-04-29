%Wrapper for AODR analysis pipeline

%Add appropriate directories
addAODRpathsADpc
AnalysisLoc = pwd;

behaviorOnly = input("Unsorted? (1/0):")
if behaviorOnly
    figLoc = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Figures\Behavior';
    monkf = input('Mr. M? (1/0):')
    if monkf
        dir('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\MrM')
    else
        dir('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Cicero')
    end
    
else
    dir('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Sorted')
    figLoc = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Figures\SortedSessions';
end

fileName = input("Paste filename:", "s")
monkey = fileName(1:2);

%Plot behavior regardless if sorted
plotADPODR_sessionBehavior(fileName, monkey, behaviorOnly, [])
f0 = gcf;
% f0.WindowState = 'maximize';
f0.Units = 'normalize';
f0.OuterPosition = [0 0 0.35 1.5];
cd(figLoc)
if ~exist([figLoc fileName],'dir')
    mkdir(fileName)
end
cd(fileName)
exportgraphics(f0,[fileName '_1_beh.png'],'Resolution',300)
close(f0)

%Plot neural data if sorted
if ~behaviorOnly
    cd(AnalysisLoc)
    neuralPlots(fileName, monkey, behaviorOnly, figLoc);
    cd(AnalysisLoc)
    neuralBehPlots(fileName, monkey, behaviorOnly, figLoc)
end