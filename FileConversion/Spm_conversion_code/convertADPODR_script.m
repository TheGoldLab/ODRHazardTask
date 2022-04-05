%% Script for converting .plx files into .mat files (FIRA) for
%   the adaptive ODR task

BoxLocal = 'C:/Users/alice'; %Alice PC

PlxLoc = '/Box/GoldLab/Data/Physiology/AODR/Sorted';
SaveLoc = '/Box/GoldLab/Data/Physiology/AODR/SortedConverted';

data_path = [BoxLocal PlxLoc];
save_path = [BoxLocal SaveLoc];

% data_path = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/Sorted';
% save_path = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/AODR/Converted';
file_in = fullfile(data_path, 'MM_2022_03_23_6_58-sort01-3units.plx'); %will try to connect to server if it can't find file, check ext
file_out = fullfile(save_path, 'MM_2022_03_23_6_58-sort01-3units');
% MM_2022_02_08_4_48-sort01-2units.plx

addpath(genpath('C:\Users\alice\Documents\MATLAB\Lab_Matlab_Utilities'));
addpath(genpath('C:\Users\alice\Documents\Projects\ODRHazardTask\FileConversion\MatlabClientDevelopKit'));
addpath(genpath('C:\Users\alice\Documents\Projects\ODRHazardTask\FileConversion\Plexon_VersionAug2013'));
addpath(genpath(pwd))

bFile(file_in, [], 'spmADPODR5', file_out, 'all', 49:51, 0, 1, 0, []);