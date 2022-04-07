%Wrapper for AODR analysis pipeline

%Add appropriate directories
addAODRpathsADpc

behaviorOnly = input("Unsorted? (1/0):")
if behaviorOnly
    dir('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\MrM')
else
    dir('C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\SortedConverted')
    figLoc = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Figures\SortedSessions';
end

fileName = input("Paste filename:", "s")
monkey = fileName(1:2);

%Plot behavior regardless if sorted
plotADPODR_sessionBehavior(fileName, monkey, behaviorOnly, axs)

if ~behaviorOnly
    neuralPlots(fileName, monkey, behaviorOnly, figLoc);
    neuralBehPlots(fileName, monkey, behaviorOnly, figLoc)
end