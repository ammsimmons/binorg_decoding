% This is an updated script that computes average decoding accuracy 
% and performes cluster-based permutation analysis. 

% In order to perform cluster-based permutation analysis, users must use
% this script to simulate null-distribution to test against. 

% This script needs all the original decoding results in the same
% directory. These simulations utilize these subject decoding results.

%Edited by Aaron Simmons (UC Davis)
%Original Author: Gi-Yeul Bae (ASU)

%% Parameters
% Must match SVM_ECOC_Stat_Plot_Permutation_2020.m 

% Subject List
subList = [505 506 507 508 509 510 512 514 516 517 519 520 521 523 524 525];       
Nsub = length(subList);

% Decoding Parameters (must match svmECOC object)
Nblock = 3; % cross-validation
Nitr = 10; % iteration
Ntp = 100; % # of time points
NBins = 16; % # of stimulus bins

% Permutation Simulation Parameters 
NPermutations = 10000;  %# of permutations
permutedT = nan(1,NPermutations); %preallocate memory
NTrials  = zeros(Nsub,1); %Preallocate trials from svmECOC object 
chancelvl = 1/NBins; %chance level of avg. decoding

% resampled timepoints from decoding (20ms/tp) 
tm = -500:20:1496;

%create empty matrix
AverageAccuracy = zeros(Nsub,Ntp);

% which time points to permute? relevantTime to conduct stats
tStart = find(tm==220); % start of retention interval
ntPerm = find(tm==1480); % end of retention interval (in this case, end of epoch)

% Load Individual subject filenames
fileLocation = pwd; 
fName = ['/Orientation_Results_ERPbased_']; 

% Save Filename (saves .mat file to current directory)
saveName = ['Permuted_T_mass.mat'];

tic

for permt = 1: NPermutations
    
for sub = 1:Nsub
    DecodingAccuracy = zeros(Ntp,Nblock,Nitr);
   
    % Load Data
    readThis =strcat(fileLocation,fName,num2str(subList(sub)),'.mat');
    load(readThis)
    
    NTrials(sub,1) = svmECOC.nTrials;
    
    svmPrediction = squeeze(svmECOC.modelPredict);
    tstTargets = squeeze(svmECOC.targets);
    
        for block = 1:Nblock
            for itr = 1:Nitr
                %% assign random target ID for Permutation Testing
                Answer = shuffle(1:NBins)'; 

                for tp = 1:Ntp % for stimulus duration of 0 ms ~ 1480 ms  
                    
                    prediction = squeeze(svmPrediction(itr,tp,block,:)); % this is predictions from models
                    Err = Answer - prediction;
                    ACC = mean(Err==0);
                    DecodingAccuracy(tp,block,itr) = ACC; % average decoding accuracy

                end
                
            end
        end
      
     grandAvg = squeeze(mean(mean(DecodingAccuracy,2),3));
     
     smoothed = nan(1,Ntp);
     for tAvg = 1:Ntp
         if tAvg ==1
           smoothed(tAvg) = mean(grandAvg((tAvg):(tAvg+2)));
         elseif tAvg ==2
           smoothed(tAvg) = mean(grandAvg((tAvg-1):(tAvg+2)));
         elseif tAvg == (Ntp-1)
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+1)));
         elseif tAvg == Ntp
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg)));
         else
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+2)));  
         end

     end
     
     
     AverageAccuracy(sub,:) =smoothed; % average across iteration and block
     
    clear svmECOC
    
end %End of subject

% AverageAccuracy = Nsub x Ncond (Attend vs Ignore) x Ntyp (Evoked VS Total) x Ntp;

%now compute average accuracy across participants
subAverage = squeeze(mean(AverageAccuracy,1)); 
seAverage = squeeze(std(AverageAccuracy,1))/sqrt(Nsub); 

%% do cluster mass analyses
 
releventTime = tStart:ntPerm; % during the retention interval 
 
Ps = nan(2,length(releventTime));
    for i = 1:length(releventTime) % make sure this time range is correct
        tp = releventTime(i);
 
        [H,P,CI,STATS] =  ttest(AverageAccuracy(:,tp), chancelvl,'tail','right'); % Test Against Zero

        Ps(1,i) = STATS.tstat;
        Ps(2,i) = P;
    end
% find significant points
candid = Ps(2,:) <= .05;

candid_marked = zeros(1,length(candid));
candid_marked(1,1) = candid(1,1);
candid_marked(1,length(candid)) = candid(1,length(candid));
%remove orphan time points
for i = 2:length(releventTime)-1
    
    if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0
    candid_marked(1,i) = 0; 
    else
    candid_marked(1,i) = candid(1,i);     
    end
    
end

% combine whole time range with relevent time & significant information
clusters = zeros(length(tm),1);
clusterT = zeros(length(tm),1);
clusters(releventTime,1) = candid_marked; % significant or not
clusterT(releventTime,1) = Ps(1,:);  % t values
clusterTsum = sum(Ps(1,logical(candid_marked)));

%%find how many clusters are there, and compute summed T of each cluster
tmp = zeros(10,25); % creates a matrix with arbitrary size (n cluster x size of each cluster)
cl = 0;
member = 0;
for i = 2:length(clusters)-1
        
        
        if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1 
        cl = cl+1;
        member = member +1;
        tmp(cl,member) = i;    
        
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0 
            if i == 2
                cl = cl +1;
                member = member +1;  
                tmp(cl,member) = i;    
                member = 0;  
                
            else
                member = member +1;  
                tmp(cl,member) = i;    
                member = 0;  
            end
        
        
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1             
            if i ==2
            cl = cl+1;
            member = member +1;  
            tmp(cl,member) = i;    
            else
            member = member +1;  
            tmp(cl,member) = i;    
            end
            
        
        else
        
            
        end
end


HowManyClusters = cl;
a = tmp(1:cl,:);
eachCluster = a(:,logical(sum(a,1)~=0));

%now, compute summed T of each cluster 
dat_clusterSumT = zeros(HowManyClusters,1);
for c = 1:HowManyClusters
   dat_clusterSumT(c,1) = sum(clusterT(eachCluster(c,eachCluster(c,:) ~=0)));
end
    
    if size(dat_clusterSumT,1) > 0 % if there is at least one signifiant cluster
   permutedT(1,permt) = max(dat_clusterSumT);
    else
   permutedT(1,permt) = 0;     
    end
if rem(permt,(NPermutations/10)) == 0   
     progress = 100*(permt/NPermutations);
     whatToPrint = strcat(num2str(progress),'% is done.');
     disp(whatToPrint)
     toc  % end of one permutation
     tic
 end

end % end of simulation

permutedT = sort(permutedT);
save(saveName,'permutedT');
