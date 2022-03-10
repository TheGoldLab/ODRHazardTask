% YF made for KS and now also for AD
data_path = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\Sorted';
save_path = 'C:\Users\alice\Box\GoldLab\Data\Physiology\AODR\SortedConverted';
convertPath = 'C:\Users\alice\Documents\Projects\ODRHazardTask\FileConversion';

% data_path = 'C:\Users\alice\Desktop\Ci_training';
% save_path = 'C:\Users\alice\Desktop\Ci_training';


% data_path = 'C:\Users\alice\Desktop\MM_training';
% save_path = 'C:\Users\alice\Desktop\MM_training';

% data_path = '/Users/lab/Desktop/MM_training/RawPLX';
% save_path = '/Users/lab/Desktop/MM_training/MatFiles';

% data_path = '/Users/lab/Desktop/Ci_training/RawPLX';
% save_path = '/Users/lab/Desktop/Ci_training/MatFiles';
% 
% convertPath = '/Users/lab/Documents/Alice/ODRHazardTask/FileConversion';
addpath(genpath(convertPath));
files = dir(fullfile(data_path, '*.plx'));
nfiles = length(files);
cd(data_path);

for i=nfiles-15%1:nfiles
    
fname_plx = files(i).name;
fname = fname_plx(1:end-4);
fname=[fname,'.mat'];
    if exist(fname)==0


bFile(fname_plx, [], 'spmADPODR4', fname, 'all',[49,50,51] , 0, 1, 0, [],save_path); %'all',{[18,49,50,51],{'LFP','pupil','eye_x','eye_y'}}
    end
  
clear fname fname_plx
end
disp('DONE!');
