function [projection cd wt] = trainCD(X,y,includeI,ops)

%% function [projection cd wt] = trainCD(X,y,includeI,ops)
%
% this function estimates the 'coding direction' between target and
% noise trials using leave-one-out crossvalidation. it does this by
% using all trials but one to compute the difference between target
% and noise trials for each cell over time (wt), then uses an
% average of these weights in a time window to determine the coding
% direction (cd). the test trial is then projected along the CD
% vector to get a population projection for each trial.
% 
% INPUTS:
%  X - PSTH with all neurons and trials (n trials x m neurons x [p timepoints])
%  y - index with the trial identity (0 = noise, >0 = target)
%  includeI - trials to include in training/testing
%  ops - ops struct
%
% ops.timeInd = logical index into time dimension of X, if X has
%               3rd dim over which CD is averaged
% ops.noiseLevel = value(s) in y to use for noise trials
% ops.targetLevel = value(s) in y to use for target trials

if isempty(includeI)
    includeI = ones(size(X,1),1);
end


% setup target and noise indices
noiseI = ismember(y,ops.noiseLevel) & includeI;
targetI = ismember(y,ops.targetLevel) & includeI;

% set iteration count
its = length(includeI);

% check for time
if size(X,3) > 1
    % if neural data has time dimension, check for time index
    if ~isfield(ops,'timeInd') | isempty(ops.timeInd)
        warning('Error in trainCD.m, no time index speficied for PSTH');
    end
else
    ops.timeInd = 1;
end
    
    

% iteration loop
for i = 1:its
    
    % all training trials, except current trial
    trainI = ones(length(includeI),1);
    trainI(i) = 0;
    
    % target vs noise trials
    targetTrainI = find(targetI & trainI);
    noiseTrainI = find(noiseI & trainI);

    % test index
    testI = ~trainI;
    
    % average over trial types
    noisePopMean = squeeze(mean(X(noiseTrainI,:,:),1));
    targetPopMean = squeeze(mean(X(targetTrainI,:,:),1));

    % weights
    wt = targetPopMean' - noisePopMean'; % imagesc(ops.psthTime,1:size(wt,1),wt);
    
    % make CD vector with the area in the target window for each
    % neuron
    CD = mean(wt(:,ops.timeInd),2);
    
    % do a projection for each test trial
    trialMat = squeeze(X(find(testI),:,:));
    projection(i,:) = trialMat*CD;
    
    % save out the CD vector to average later
    cd(i,:) = CD;
    
end



