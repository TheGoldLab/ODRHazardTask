%Get the file
data_path = '/Users/lab/Desktop/MM_training/MatFiles';
figPath = '/Users/lab/Desktop/MM_training/Figures';
files = dir(fullfile(data_path, '*.mat'));
nfiles = length(files);
cd(data_path);

currentFile = 7;

theFile= files(currentFile).name;

cd(figPath)
[~,b,~]=fileparts(theFile);
mkdir(b)