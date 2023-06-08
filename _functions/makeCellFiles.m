clear all; close all;

warning off

%% directory stuff
behaviorBase = '~/gits/gain-gonogo';
spikeBase = '~/data/kilosort';
wd = '~/chris-lab/projects/contrast_behavior';
addpath(genpath(behaviorBase));
addpath(genpath('~/gits/npy-matlab'));
addpath(genpath('~/chris-lab/code_general'));
addpath(genpath([wd filesep]));

% mouse and saving options
mouseList = {'CA051','CA061','CA070','CA073','CA074',...
             'CA102','CA104','CA105','CA106','CA107',...
             'CA118','CA119','CA121','CA122'};
datDir = '~/data/gain_behavior/_sessions';
stimDir = '~/data/gain_behavior/_stimuli';

% match neural to behavioral data
spikesess = matchBehaviorToSpikes(mouseList,behaviorBase,spikeBase,30);

% files to skip
skipList(1,:) = {29,'empty event ttls'};
skipList(2,:) = {68,'events dont align'};
skipList(3,:) = {173,'events stop early'};

%% session loop
cnt = 0;
cells = [];
spikes = [];
cellinfo = [];
clear session behavior events waveform;
for i = 1:length(spikesess)

    % if this isn't a skip session
    if sum(i == [skipList{:,1}]) == 0
        
        
        % preprocess the data file
        fprintf('MOUSE %s - SESSION %s\n',spikesess(i).mouse,spikesess(i).session);
        fnout = makeSessionData(spikesess(i).mouse,spikesess(i), ...
                                datDir);
        
        if ~isempty(fnout)
            
            cnt = cnt + 1;
            
            % load the data
            d = load(fnout);
                        
            %% cell data
            spikes = [spikes; d.spikes'];
            cellinfo = [cellinfo; d.data.cellInfo];
            
            % waveform info:
            fields2remove = {'chanMap','waveform','max','noise_level',...
                             'percentile95','snr','stddev','peak','fs'};
            if cnt == 1
                waveform = rmfield(d.data.wave,fields2remove);
            else
                fields = fieldnames(waveform);
                pruned_waveform = rmfield(d.data.wave,fields2remove);
                for j = 1:length(fields)
                    waveform.(fields{j}) = [waveform.(fields{j}); ...
                                        pruned_waveform.(fields{j})];
                end
            end
            disp(waveform)
            
            %% session data
            session(cnt) = d.sessionInfo;
            behavior(cnt) = d.data.behavior;
            events(cnt) = d.data.events;
            
        end

    end
    
end

% save spike data
save(fullfile(datDir,'..','spikeData.mat'),'spikes','cellinfo','waveform','-v7.3');

% save session data
save(fullfile(datDir,'..','sessionData.mat'),'session','behavior','events');

% for each mouse, visualize number of neurons over time
ncells = vertcat(session.nCells);
dates = vertcat(session.sessionID);
mice = cellstr(vertcat(session.mouse));

figure;
% get dates
dt = datetime(num2str(dates),'InputFormat','yyMMddHHmm')
uM = unique(mice);
for i = 1:length(uM)
    
    I = contains(mice,uM{i});
    dvec = round(days(dt(I) - min(dt(I))));
    ncvec = ncells(I);
    
    hold on;
    plot(dvec,ncvec);
    
end
