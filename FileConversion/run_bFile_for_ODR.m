% YF made for KS and now also for AD
data_path = '/Users/lab/Desktop/MM_training/RawPLX';
files = dir(fullfile(data_path, '*.plx'));
nfiles = length(files);
cd(data_path);

for i=1:nfiles
    
fname_plx = files(i).name;
fname = fname_plx(1:end-4);
fname=[fname,'.mat'];
    if exist(fname)==0


bFile(fname_plx, [], 'spmADPODR3', fname, 'all', [49,50,51], 0, 1, 0, []);
end
   
end
disp('DONE!');
