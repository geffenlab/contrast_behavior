clear all
close all

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


session = matchBehaviorToSpikes(mouseList,behaviorBase,spikeBase,30);

skipList(1,:) = {29,'empty event ttls'};
skipList(2,:) = {68,'events dont align'};
skipList(3,:) = {173,'events stop early'};

cnt = 0;
cells = [];
for i = 1:length(session)

    if (sum(i == [skipList{:,1}]) == 0)
        
        fprintf('MOUSE %s - SESSION %s\n',session(i).mouse,session(i).session);
        fnout = makeSessionData(session(i).mouse,session(i),fullfile(datDir));
        
        % load the data
        if ~isempty(fnout)
            load(fnout);
        end
        
        %b = load(data.behavior.matFile);
        %if isfield(b,'level')
        %    keyboard
        %end
        
        % check for stimulus data
        fp = strsplit(data.params.stimInfo.stim,'\');
        stimFile = fullfile(stimDir,fp{end});
        
        % NOTE: on Jan 29, 2019 I found a bug where the level index
        % on noise trials was unknown... need to test to see if I
        % can recover data from earlier files...
        % https://github.com/geffenlab/gain-gonogo/commit/5870159e1e0df2d92c7b2075391844733c40faf2
        
        if exist(stimFile,'file') & isfield(data.behavior,'level') & sum(data.behavior.level) > 0
                        
            % load the stim file
            stim = load(stimFile,'DB','params');
            
            if isfield(stim,'DB')
                % analysis of each cell
                for c = 1:size(spikes,2)
                    
                    % mean FR
                    fr = length(spikes{c}) / ...
                         (data.events.trialOn(end) - data.events.trialOn(1));
                    
                    % use only cells with FR greater than 1Hz
                    if fr > 1
                        
                        cnt = cnt + 1;
                        
                        % multi or single unit -- ISI distribution
                        contamPct = sum(diff(spikes{c}) < .002) ./ ...
                            length(spikes{c});
                        if round(contamPct,2) > .01
                            lbl = 'mu';
                        else
                            lbl = 'su';
                        end
                        
                        % PSTH at chord resolution
                        edges = [-.5:.005:5.5];
                        ts = edges(2:end) - mean(diff(edges))/2;
                        trialEvs = data.events.trialOn;
                        targEvs = data.events.targOn;
                        psth = makePSTH(spikes{c},trialEvs,edges);
                        
                        [~,sortI] = sortrows(data.behavior.trialType(:,3));
                        
                        % compute the noise ratio for each contrast only during
                        % noise trials, excluding responses near background changes
                        noiseI = data.behavior.trialType(:,1)<1;
                        npsth = psth(noiseI,:);
                        ind(1,:) = ts > .3 & ts < 3;
                        ind(2,:) = ts > 3.3 & ts < 4;
                        NR = [];
                        for i = 1:2
                            NR(:,i) = responsePower(npsth(:,ind(i,:)),data.behavior.trialType(noiseI,3))';
                        end
                        
                        % compute waveform info
                        if ~isempty(data.wave.waveform)
                                                        
                            wave.chanMap = data.wave.chanMap;
                            wave.chan_waveforms = data.wave.waveform;
                            wave.peak_chan = data.wave.peakchan(c);
                            wave.waveform = data.wave.peakwaveform(c,:);
                            wave.spikecount = data.wave.spkcount(c);
                            wave.snr = data.wave.snr(c,wave.peak_chan);
                            wave.trough_peak = data.wave.trough_peak(c);
                            wave.peak_inflect = data.wave.peak_inflect(c);
                            wave.FWHM = data.wave.FWHM(c);
                            
                        else
                            
                            wave.chanMap = nan;
                            wave.chan_waveforms = nan;
                            wave.peak_chan = nan;
                            wave.waveform = nan(1,63);
                            wave.spikecount = nan;
                            wave.snr = nan;
                            wave.trough_peak = nan;
                            wave.peak_inflect = nan;
                            wave.FWHM = nan;
                            
                        end
                        
                        
                        
                        
                        % add in data for each neuron
                        cells(cnt).mouse = data.params.mouse;
                        cells(cnt).session = num2str(data.params.session.sessionInfo(end));
                        cells(cnt).cellnum = data.cellInfo{c,4};
                        cells(cnt).cellcnt = c;
                        cells(cnt).cellType = lbl;
                        cells(cnt).cellID = sprintf('%s-%s-%02d-%03d-%s',...
                                                    cells(cnt).mouse,...
                                                    cells(cnt).session,...
                                                    cells(cnt).cellcnt,...
                                                    cells(cnt).cellnum,...
                                                    cells(cnt).cellType);
                        cells(cnt).spikes = spikes{c};
                        cells(cnt).trialEvents = trialEvs;
                        cells(cnt).targetEvents = targEvs;
                        cells(cnt).lickEvents = data.events.lick;
                        cells(cnt).rewardEvents = data.events.rewd;
                        %cells(cnt).psth = psth;
                        %cells(cnt).psthTS = ts;
                        cells(cnt).NR = NR;
                        cells(cnt).chanMap = wave.chanMap;
                        cells(cnt).chan_waveforms = wave.chan_waveforms;
                        cells(cnt).peak_chan = wave.peak_chan;
                        cells(cnt).waveform = wave.waveform;
                        cells(cnt).spikecount = wave.spikecount;
                        cells(cnt).snr = wave.snr;
                        cells(cnt).trough_peak = wave.trough_peak;
                        cells(cnt).peak_inflect = wave.peak_inflect;
                        cells(cnt).FWHM = wave.FWHM;
                        cells(cnt).ISIcontamination = contamPct;
                        cells(cnt).stimFile = stimFile;
                        cells(cnt).behaviorFile = data.params.session.behavior;
                        cells(cnt).neuralFile = fnout;
                        
                        % add in data about the stimulus/behavior
                        %cells(cnt).spec = stim.DB;
                        %cells(cnt).stimInfo = stim.params;
                        cells(cnt).offsets = data.params.offsets;
                        cells(cnt).trialType = data.behavior.trialType;
                        cells(cnt).response = data.behavior.response;
                        cells(cnt).abort = data.behavior.abort;  
                        cells(cnt).level = data.behavior.level;
                        
                        % mean noise ratio
                        nr = cells(cnt).NR;
                        nr(nr<0) = nan;
                        nr = nanmedian(nr,1);
                        
                        % save the file
                        spec = stim.DB;
                        stimInfo = stim.params;
                        d = cells(cnt);
                        fn = sprintf('%s.mat',cells(cnt).cellID);
                        fp = fullfile('~/data/gain_behavior','_cells',fn);
                        fprintf('\tSaving cell %s\n',fn);
                        save(fp,'d','spec','stimInfo');
                        
                    end
                end
            else
                fprintf(['\Stimfile doesnt have DB field... ' ...
                         'skipping\n']);
            end
        else
            fprintf('\tNo spectrogram found... skipping\n');
        end
    else
        fprintf('\tFile in skip list... skipping\n');
    end
end

keyboard

%cells = rmfield(cells,'psth');
%cells = rmfield(cells,'psthTS');

save(fullfile(datDir,'contrastGLM','_allCellsGLM.mat'),'cells','-v7.3');

fields = fieldnames(cells(1));
for i = 1:length(fields)
    tmp = {cells.(fields{i})};
    pp = whos('tmp');
    fprintf('%3.2fMB - %s\n',pp.bytes*1e-6,fields{i});
end


%% prune data by noise ratio
clear NR;
for i = 1:length(cells)
    nr = cells(i).NR;
    nr(nr<0) = nan;
    NR(i,:) = nanmedian(nr,1);
end

% only include cells with low noise ratios in BOTH contrasts
NRcut = 50;
nrI = all(NR < NRcut,2);
nrSub = find(nrI);

fields = fieldnames(cells(1));

neuron = [];
for i = 1:length(nrSub)
    for j = 1:length(fields)
        neuron(i).(fields{j}) = cells(nrSub(i)).(fields{j});
    end
end

% find unique stim files for these neurons
stimFiles = unique({neuron.stimFile});

% save an aggregate file
save(fullfile(datDir,'contrastGLM','_allCellsGLM_lowNR.mat'), ...
     'neuron','stimFiles','-v7.3');

keyboard





