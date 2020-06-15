% This is an updated function of Bae & Luck (2018) Decoding
% pipeline that utilizes a nested bin-epoched data structure.
% Refer to OSF: https://osf.io/29wre/

% NOTE: low-pass filtering to 6hz was applied to continuous EEG (prior to
% binning)


% Edited by Aaron Simmons


function SVM_ECOC_ERP_Decoding(subs)
% delete(gcp)
% parpool

if nargin ==0
    
    subs = [505];
    
end

nSubs = length(subs);


% parameters to set

svmECOC.nChans = 16; % # of channels

svmECOC.nBins = svmECOC.nChans; % # of stimulus bins

svmECOC.nIter = 10; % # of iterations

svmECOC.nBlocks = 3; % # of blocks for cross-validation

svmECOC.frequencies = [0 6]; % low pass filter

svmECOC.dataTime.pre = -500; % Set epoch start (from imported data)
svmECOC.dataTime.post = 1496; % Set epoch end (from imported data)

svmECOC.time = -500:20:1496; % time points of interest in the analysis

svmECOC.window = 4; % 1 data point per 4 ms in the preprocessed data

svmECOC.Fs = 250; % samplring rate of in the preprocessed data for filtering

ReleventChan = sort([2,3,4,18,19, 5,6,20, 7,8,21, 9,10,11,12,13,14, 22,23,24,25,26,27, 15,16,1,17]); %electrodes

svmECOC.nElectrodes = length(ReleventChan); % # of electrode included in the analysis


% for brevity in analysis

nChans = svmECOC.nChans;

nBins = svmECOC.nBins;

nIter = svmECOC.nIter;

nBlocks = svmECOC.nBlocks;

freqs = svmECOC.frequencies;

dataTime = svmECOC.dataTime;

times = svmECOC.time;

nElectrodes = svmECOC.nElectrodes;

nSamps = length(svmECOC.time);

Fs = svmECOC.Fs;

%% Loop through participants
for s = 1:nSubs
    sn = subs(s);
    
    fprintf('Subject:\t%d\n',sn)
    
    % load data
    currentSub = num2str(sn);
    dataLocation = pwd; % set directory of data set
    loadThis = strcat(dataLocation,'/Decoding_',currentSub,'_BE.mat');
    load(loadThis)
    
    % where to save decoding output
    saveLocation = pwd; % set directory for decoding results.
    
    % grab EEG data from bin-list organized data 
    eegs = binorgEEG.binwise_data;
    nPerBin = binorgEEG.n_trials_this_bin;
    
    % set up time points
    tois = ismember(dataTime.pre:4:dataTime.post,svmECOC.time); nTimes = length(tois);
    
    % # of trials
    svmECOC.nTrials = (sum(nPerBin)); nTrials = svmECOC.nTrials;
    
    
    % Preallocate Matrices
    
    svm_predict = nan(nIter,nSamps,nBlocks,nChans); % a matrix to save prediction from SVM
    tst_target = nan(nIter,nSamps,nBlocks,nChans);  % a matrix to save true target values
    
    % create svmECOC.block structure to save block assignments
    svmECOC.blocks=struct();
    
    if nChans ~= max(size(eegs))
        error('Error. nChans dne # of bins in dataset!');
        % Code to be used on a bin-epoched dataset for one
        % class only (e.g. Orientation). Users may
        % try to bin multiple classes in one BDF.
        % This is not currently allowed.
    end
    
    
    % Loop through each iteration
    tic % start timing iteration loop
    
    for iter = 1:nIter

        %using bin-epoched datasets 
        %preallocate      
        blockDat_filtData = nan(nBins*nBlocks,nElectrodes,nSamps);   
        labels = nan(nBins*nBlocks,1);     % bin labels for averaged & filtered EEG data
        blockNum = nan(nBins*nBlocks,1);   % block numbers for averaged & filtered EEG data
        bCnt = 1; %initialize binblock counter 
        
        for bin = 1:nChans
            
            % We will use as many possible trials per
            % bin having accounted already for artifacts
            
            %Drop excess trials 
            nPerBlock = floor(nPerBin/nBlocks); %array for nPerBin
            
            %Obtain index within each shuffled bin
            shuffBin = randperm((nPerBlock(bin))*nBlocks)';
            
            %Preallocate arrays
            
            blocks = nan(size(shuffBin));
            shuffBlocks = nan(size(shuffBin));
            
            %arrage block order within bins
            x = repmat((1:nBlocks)',nPerBlock(bin),1);
            shuffBlocks(shuffBin) = x;
            
            %unshuffled block assignment
            blocks(shuffBin) = shuffBlocks;
            
            % save block assignment
            blockID = ['iter' num2str(iter) 'bin' num2str(bin)];
            
            svmECOC.blocks.(blockID) = blocks; % block assignment
            
            
            %create ERP average and place into blockDat_filtData struct          
            %grab current bin with correct # of electrodes & samps
            eeg_now = eegs(bin).data(ReleventChan,:,:);
            
            for bl = 1:nBlocks
                
                blockDat_filtData(bCnt,:,:) = squeeze(mean(eeg_now(:,tois,blocks==bl),3));
                
                labels(bCnt) = bin;
                
                blockNum(bCnt) = bl;
                
                bCnt = bCnt+1;
                
            end
            
        end
        
        
        % Do SVM_ECOC at each time point
        parfor t = 1:nSamps
            
            % grab data for timepoint t
            
            toi = ismember(times,times(t)-svmECOC.window/2:times(t)+svmECOC.window/2);
            
            % average across time window of interest
            
            dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));
            
            % Do SVM_ECOC for each block
            for i=1:nBlocks % loop through blocks, holding each out as the test set
                
                
                trnl = labels(blockNum~=i); % training labels
                
                tstl = labels(blockNum==i); % test labels
                
                trnD = dataAtTimeT(blockNum~=i,:);    % training data
                
                tstD = dataAtTimeT(blockNum==i,:);    % test data
                
                % SVM_ECOC
                mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM' );   %train support vector mahcine
                
                LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                
                svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                
                tst_target(iter,t,i,:) = tstl;             % save true target labels
                
                
            end % end of block
            
        end % end of time points
        
    end % end of iteration
    
    toc % stop timing the iteration loop
    
    
    
    OutputfName = strcat(saveLocation,'/Orientation_Results_ERPbased_',currentSub,'.mat');
    
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict;
    
    svmECOC.nBlocks = nBlocks;
    
    save(OutputfName,'svmECOC','-v7.3');
    
end % end of subject loop