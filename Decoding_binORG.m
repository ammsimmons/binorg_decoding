% binorgEEG script:
% Use *right after* Extracting Bin-Based Epochs
% Creates bin-organized .mat files for Decoding
% Writes .mat files to current working directory

% Place continuous BE datasets (.set & .fdt) in current working directory
% Examples of these kind of data can be found on box: 
% https://ucdavis.box.com/s/jrd9zus448nkc0wjkjf1x5nxbc7oxwg8

% by Aaron M. Simmons, University of California, Davis

clear all;
close all; 


parentfolder = pwd; 
subject_list = {505, 506, 507, 508, 509, 510, 512, 514, 516, 517, 519, 520, 521, 523, 524, 525};
numsubjects = length(subject_list);


for s = 1:numsubjects
    subject = subject_list{s};
    subjectfolder = [parentfolder]; %loc of file
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, num2str(subject));
    
    %Initialize EEG
    eeglab;
    
    %load data
    EEG = pop_loadset('filename', [num2str(subject) '_binned_be.set'], 'filepath', subjectfolder);
    
    %use binorg-EEG function (this function is only on ERPlab V8.1)
    binepEEG_to_binorgEEG(EEG, ['Decoding_BE_' num2str(subject)]); 
    % Produces bin-organized .mat file and outputs file into current working directory
    % it will be called "Decoding_BE_XXX" as specified in string filename
   
    close all;
    clear EEG;
    
end

