% binorgEEG script:
% Use *right after* Extracting Bin-Based Epochs
% Creates bin-organized .mat files for Decoding
% Writes .mat files to current working directory

% Place original EEGlab datasets in current working directory

% by Aaron M. Simmons, University of California, Davis

clear all;
close all; 


parentfolder = pwd; 
subject_list = {505, 506, 507, 508, 509, 510, 512, 514, 516, 517, 519, 520, 523, 524, 525};
numsubjects = length(subject_list);


for s = 1:numsubjects
    subject = subject_list{s};
    subjectfolder = [parentfolder]; %loc of file
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, num2str(subject));
    
    %Initialize EEG
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %load data
    EEG = pop_loadset('filename', [num2str(subject) '_binned_be.set'], 'filepath', subjectfolder);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [num2str(subject)], 'gui', 'off');
    
    %use binorg-EEG function (this function is only on ERPlab V8.1)
    binepEEG_to_binorgEEG(EEG, ['Decoding_BE_' num2str(subject)]); 
    %writes to current directory 
   
    close all;
    clear eeglab;
    
end

eeglab redraw;
erplab redraw;
