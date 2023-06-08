function [stim,spec,y0,ops] = formatStimAndSpikes(spikes,events,stimInfo,behavior,ops,flag)

% file
specFile = fullfile('./_data','_spectrograms',[stimInfo.IDsess '-spec.mat']);

% spectrogram file isnt found, make it
if ~exist(specFile,'file')
    
    % load spectrograms
    fp = strsplit(stimInfo.stim,'\');
    fprintf('\tLoading stimuli... '); tic;
    S = load(fullfile(pwd,'_stimuli',fp{end})); toc;
    
    % parse the stimulus
    for i = 1:length(behavior.trialType)
        % index the stimulus structure for this trial
        tt = behavior.trialType(i,:);
        if tt(1) > 0
            ind = tt(1);
        else
            ind = 1;
        end
        
        stimIndex(i,:) = [(tt(1)>0)+1,tt(2),tt(3),ind];
        
        % stim length in model samples, with just the stimulus
        stimLength(i) = length(S.DB{stimIndex(i,1),stimIndex(i,2),stimIndex(i,3),stimIndex(i,4)}) * ...
            (ops.fs * stimInfo.chordDuration);
        
        % stim length plus pads
        padStimLength(i) = 2*(ops.pad * ops.fs) + stimLength(i);

    end
    stim_samps = sum(stimLength);
    nsamps = sum(padStimLength);
    minStimLen = min(stimLength);
    minPadStimLen = min(padStimLength);
    
    % pad samples
    ps = ops.pad * ops.fs;
    
    % preallocate the full stim
    stim.db = zeros(length(ops.f),nsamps);
    stim.index = zeros(7,nsamps);
    ind = 0;
    
    % for each trial, build the stimulus and index
    for i = 1:length(behavior.trialType)
        
        % trial stimulus spectrogram
        St = cleanResample(S.DB{stimIndex(i,1),stimIndex(i,2),...
                            stimIndex(i,3),stimIndex(i,4)},...
                           stimInfo.chordDuration,1/ops.fs);
        
        % number of samples for this stimulus
        ns = min([ops.trialCutoff length(St)]);
        
        % number of samples including pad
        n = ns + 2 * ps;
        
        % trial index
        ind = ind(end) + 1:ind(end)+n;
        
        % target index
        targT = stimInfo.offsets(behavior.trialType(i,2)) + ...
                stimInfo.baseNoiseD;
        target = targT * ops.fs;
        
        
        
        %% stimulus (with padding on both ends)
        stim.db(:,ind) = [zeros(length(ops.f),ps) St(:,1:ns) zeros(length(ops.f),ps)];
        
    
        clear index;
        %% trial index
        index(1,:) = ones(1,ns)*i;
        
        
        %% contrast index
        index(2,:) = [ones(1,stimInfo.baseNoiseD*ops.fs) ...
                    ones(1,ns-stimInfo.baseNoiseD*ops.fs)*2];
        
        
        %% stimulus pattern
        index(3,:) = ones(1,ns) * behavior.trialType(i,3);
        
        
        %% switches (trial start and contrast switch)
        index(4,:) = [zeros(1,ops.w*ops.fs) ...
                    ones(1,(stimInfo.baseNoiseD-ops.w)*ops.fs) ...
                    zeros(1,ops.w*ops.fs)...
                    ones(1,ns - (stimInfo.baseNoiseD+ops.w)*ops.fs)];

        
        %% after the target
        index(5,:) = [ones(1,target) ...
                      ones(1,(ops.w+ops.targetExclude)*ops.fs) .* behavior.trialType(i,1) == 0 ...
                      ones(1,ns - (target + (ops.w+ops.targetExclude)*ops.fs))];

        
        %% minimum stimulus length
        index(6,:) = [ones(1,minStimLen) zeros(1,ns-minStimLen)];
 
        
        %% compute matching windows for each contrast
        % find samples after the switch that are not after switch
        % or target
        ts = index(5,:) & index(4,:) & index(2,:)==2;
        
        % set samples immediately before the switch to be included
        pres = [zeros(1,stimInfo.baseNoiseD*ops.fs - sum(ts)) ...
                ones(1,sum(ts)) zeros(1,ns - stimInfo.baseNoiseD*ops.fs)];
        
        % add to the index
        index(7,:) = ts | pres;
        
        % check
        if sum(find(index(7,:))<=(3*ops.fs)) ~= sum(find(index(7,:))>3*ops.fs)
            fprintf('WINDOW MISMATCH');
            keyboard
        end
        
        
        %% pad all the indices (trial index extends into the pad)
        pad = [ones(1,ps)*i; zeros(size(index,1)-1,ps)];
        stim.index(:,ind) = [pad index pad];
        
        % actual stim length (cutoff)
        stimLengthUsed(i) = ns;
        
    end
    
    % trim down the indices based on cut off stimuli used
    actual_nsamps = sum(stimLengthUsed) + 2 * ps * length(behavior.trialType);
    stim.db = stim.db(:,1:actual_nsamps);
    stim.index = stim.index(:,1:actual_nsamps);
    
    % save out the data
    spec.trial_actual_length = stimLengthUsed;
    spec.trial_index = stimIndex;
    spec.trial_length = stimLength;
    spec.pad_trial_length = padStimLength;
    spec.nsamps = nsamps;
    spec.total_samps = stim_samps;
    spec.total_padded_samps = nsamps;
    spec.total_actual_samps = actual_nsamps;
    spec.min_trial_length = minStimLen;
    spec.min_padded_trial_length = minPadStimLen;
    spec.pad_samps = ps;
    spec.spec_file = specFile;
    
    save(specFile,'stim','spec');
    
else
    load(specFile);
        
end

% make spike count histogram over the stimulus
fprintf('\tFormatting spikes... '); tic;
y0 = zeros(1,length(stim.db));

ind = 0;
for i = 1:length(events)
    
    % index
    n = spec.trial_actual_length(i) + 2*spec.pad_samps;
    ind = ind(end) + 1:ind(end)+n;
    
    % get spikes ops.pad before and after the stimulus for this trial
    edges = -ops.pad:stimInfo.chordDuration:((n-spec.pad_samps)/ops.fs);
    spks = spikes - events(i);
    y0(ind) = histcounts(spks,edges);
    
end
toc;


