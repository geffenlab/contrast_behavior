function [fnout] = makeSessionData(mouse,session,output)

%% function fnout = makeSessionData(mouse,session,output)
%
% This function formats data from Chris Angeloni's gain control
% behavioral go-nogo tasks with simultaneously recorded
% neurons.
%
% It saves a mat file in the directory specified by output, with
% the following contents:
%  - trial: a struct where each entry is a trial and each 
%           subfield contains information about that trial. In 
%           most cases, the information is a timestamp relative 
%           to the start of the trial, but there are also categorical
%           variables detailing the trial type, behavioral choice, etc.
%  - spikes: a struct where each element contains spike times for a neuron
%  - cellInfo: a struct where each row is a cell, and each column 
%              contains information about that cell
%  - event: a struct with a field for each event type
%  - behavior: a struct with information about the behavioral task and
%              each trial
%  - params: a struct with general recording information

% check if output file exists
%  taskTypes = {'training','psychometric','offset'};
%  task = taskTypes{session.sessionInfo(2)};
task = session.type;
fnout = fullfile(output,sprintf('%s-%s-%s.mat',mouse, ...
                                num2str(session.sessionInfo(6)), ...
                                task));
data = [];

if ~exist(fnout,'file')
    
    fprintf('MOUSE %s - %s\n',mouse, ...
            num2str(session.sessionInfo(6)));
    

    %%  setup directories
    % spike data root directory
    root = fullfile(session.neural,...
                    'experiment1','recording1');
    
    % for openEphys format
    if ~exist(root,'dir')
        root = [root filesep '..' filesep '..'];
    end
    [broot,bf,ext] = fileparts(session.behavior);
    
    
    
    %% if there are non noise units, cancel analysis
    % check for cluster groups
    nodataflag = true;
    folders = dir(fullfile(root,'continuous'));
    if any(contains({folders.name},'Rhythm_FPGA-100.0'))
        recFolder = 'Rhythm_FPGA-100.0';
    elseif any(contains({folders.name},'Bandpass_Filter-107_100.0'))
        recFolder = 'Bandpass_Filter-107_100.0';
    else
        nodataflag = false;
    end
    
    if nodataflag
        d = fullfile(root,'continuous',recFolder);
        if exist([d filesep 'cluster_groups.csv'],'file')
            fn = [d filesep 'cluster_groups.csv'];
            tmp = importdata(fn,',');
            labels = vertcat(tmp(2:end));
        elseif exist([d filesep 'cluster_group.tsv'],'file')
            fn = [d filesep 'cluster_group.tsv'];
            tmp = importdata(fn,'\t');
            labels = vertcat(tmp(2:end));
        else
            nodataflag = false;
            labels = 'noise';
        end
    else
        labels = 'noise';
    end

    
    % continue only of there are some non-noise units
    if ~all(contains(labels,'noise'))

        
        %% BEHAVIOR INFO
        bMat = session.behavior;
        bTxt = dir([broot filesep bf '.txt']);
        bTxt = bTxt([bTxt.bytes] > 1000);

        % load it for regular trial info and determine task type
        b = load(bMat);
        
        % add abort vector if none exists
        if ~isfield(b,'abort')
            b.abort = zeros(length(b.resp),1);
        end
        if ~isfield(b,'level')
            b.level = zeros(length(b.resp),1);
        end
        
        % take miminum of behavior vectors
        [mn,mi] = min([length(b.tt) length(b.resp) length(b.abort) length(b.level)]);
        response = b.resp(1:mn);
        trialType = b.tt(1:mn,:);
        abort = b.abort(1:mn);
        level = b.level(1:mn);
        
        
        % reshape to columns
        if size(abort,1) == 1
            abort = abort';
        end
        if size(response,1) == 1
            response = response';
        end
        
        % find trials where the mouse licked
        [~,~,~,~,goodTrials] = computePerformanceGoNoGo(response,...
                                                        trialType(:,1), ...
                                                        [],[],.1);
        

        % get the time of each lick relative to the trial onset
        % (first, need to parse the text file)
        fn = [bTxt(end).folder filesep bTxt(end).name];
        [t,~,~,~,lickTimes,ts,tags] = parseLogNeural(fn);
        for i = 1:length(t)
            lickTime{i} = (t(i).lickTimes - t(i).stimTimes) / 1e6;
            if ~isempty(t(i).rewardTimes)
                rewdTime{i} = (t(i).rewardTimes(1) - t(i).stimTimes) / 1e6;
            else
                rewdTime{i} = [];
            end
        end
        
        % write out to a trial struct
        for i = 1:length(t)
            % for each trial, assign behavioral times
            trial(i).lickTime = lickTime{i};
            trial(i).rewdTime = rewdTime{i};
            trial(i).response = response(i);
            trial(i).abort = abort(i);
            trial(i).trialType = trialType(i,:);
            trial(i).goodTrial = goodTrials(i);
        end





        %% EVENT STUFF
        % this is a shitshow.
        recInfoFile = fullfile(root,'recInfo.mat');
        [ev, blockStart, startTime, msgtext,fs] = ...
            loadEventData(root,recInfoFile);

        % get the start time of this behavioral block in secs and remove
        % all events before
        stimOn = ev.ts(ev.state == 1 & ev.ts >= blockStart);
        stimOff = ev.ts(ev.state == -1 & ev.ts >= blockStart);
        lickOn = ev.ts(ev.state == 2 & ev.ts >= blockStart);
        lickOff = ev.ts(ev.state == -2 & ev.ts >= blockStart);
        rewdON = ev.ts(ev.state == 3 & ev.ts >= blockStart);
        rewdOFF = ev.ts(ev.state == -3 & ev.ts >= blockStart);
        
        % extract neural and behavioral timestamps relative to the
        % start
        respwin = [t.respWin];
        respwin = respwin(1,:)';
        bts = sort([[t.stimTimes]'; respwin])/1e6;
        bts = bts - bts(1);
        nts = stimOn - stimOn(1);
        trialID = repmat([1:length(t)],2,1);
        trialID = trialID(:);
        
        
        % plot the neural and behavior events we're starting with
        f123 = figure(123); clf; 
        subplot(1,5,[1 2]); hold on;
        plot(bts(1:end-1),diff(bts),'k')
        plot(nts(1:end-1),diff(nts),'gx')
        legend('Arduino events','recording events');
        xlabel('event time'); ylabel('event diff');
        title(sprintf('INITIAL'));
        
        
        
        % if there are large event differences early, drift correct first
        if any(diff(bts(1:10))>100)
            scales = linspace(.9,1.1,10000);
            for i = 1:length(scales)

                % mse of central chunk
                chunk = 30:30+length(bts)/2;
                mserr(i) = immse(bts(chunk)*scales(i),nts(chunk));
                
            end     
            [minerr,mI] = min(mserr);
            
            % rescale the behavior
            bts = bts * scales(mI);
            bts = bts - bts(1);
        end
        
        
        
        %% get rid of odd events
        % remove short reward events
        rewdOn = rewdON((rewdOFF-rewdON) > .005);
        rewdOff = rewdOFF((rewdOFF-rewdON) > .005);
        
        % remove events that are too short in duration (maybe they
        % were cut short by a time out)
        if ~isempty(stimOff)
            if any((stimOff - stimOn)<.0049)
                fprintf('Short event detected, removing... \n');
                ind = find((stimOff - stimOn)<.0049);
                stimOn(ind) = [];
                stimOff(ind) = [];
            end
        end
        
        
        
        
        %%  simpler removal of extra events early in the session
        % 1. algorithm looks for the first large time difference
        %    between behavior and neural events in the first nT trials
        % 2. removes that event and repeats
        
        nT = 8;
        cnt = 1;
        errI = find(abs(bts(1:nT)-nts(1:nT))>.1,1,'first');
        while ~isempty(errI)
            stimOn(errI) = [];
            stimOff(errI) = [];
            nts = stimOn - stimOn(1);
            errI = find(abs(bts(1:nT)-nts(1:nT))>.1,1,'first');
            if cnt > 3
                fprintf(['makeSessionData.m line 231: stuck removing ' ...
                         'early events!!']);
                keyboard
            end
            cnt = cnt + 1;
        end
        
       
        
        %% remove trial onsets whose differences are impossibly short
        % (possibly a buffering error? in either openEphys or
        % Matlab... quite rare in any case...)
        ITIcut = 1.5;
        if any(diff(stimOn)<ITIcut)
            fprintf('Short inter-trial-interval detected, removing...\n');
            ind = find(diff(stimOn)<ITIcut)+1;
            stimOn(ind) = [];
            stimOff(ind) = [];
        end
        nts = stimOn - stimOn(1);
                
                
        
        %% relative drift correction and event matching
        %
        % this section does drift correction and event matching:
        %  1) Correct for possible arduino clock reset
        %  2) Brute force, looks for best linear rescaling of the
        %     behavior events that matches the neural events, then
        %     finds potential error events and removes them
        %  3) Fixes the output event structure
        
        % 1)
        % check for arduino clock reset
        clockReset = find(diff(bts)>1000,1,'first');
        if ~isempty(clockReset)
            % use neural events to estimate true event diff
            truedt = nts(clockReset+1) - nts(clockReset);
            
            % reset difference
            resetdt = bts(clockReset+1) - bts(clockReset) - truedt;
            
            % correct the rest of the events
            bts(clockReset+1:end) = bts(clockReset+1:end) - resetdt;
            
        end
             
        % 2)
        % rescale and adjust bad events
        errI = [];
        errev1 = 1;
        cnt = 1;
        while ~isempty(errev1)
            
            % for different rescalings of behavior events, find the lowest mse
            scales = linspace(.9,1.1,100000);
            for i = 1:length(scales)

                % mse of central chunk (removing long events)
                chunk = 30:30+length(bts)/2;
                mserr(i) = immse(bts(chunk)*scales(i),nts(chunk));
                
            end     
            [minerr,mI] = min(mserr);
            
            % rescale the behavior
            bts = bts * scales(mI);
            bts = bts - bts(1);
            
            % plot the error
            f123 = figure(123); clf; 
            subplot(1,5,3); hold on;
            plot(scales,log(mserr),'k');
            plot(scales(mI),log(mserr(mI)),'rx');
            text(.1,.1,sprintf('drift_m = %7.8f',scales(mI)),...
                 'Units','normalized');
            xlabel('drift factor'); ylabel('log(mse)'); axis tight;
            
            % plot the neural and rescaled behavior events
            subplot(1,5,[1 2]); hold on;
            plot(bts(1:end-1),diff(bts),'k')
            plot(nts(1:end-1),diff(nts),'gx')
            legend('Arduino events','recording events');
            xlabel('event time'); ylabel('event diff');
            title(sprintf('PASS %d',cnt));

            % plot events that are quite different
            mn1 = min([length(bts),length(nts)]);
            bev = bts(1:mn1);
            nev = nts(1:mn1);
            tID = trialID(1:mn1);
            truedev = bev-nev;
            dev = abs(truedev);
            errev = nev(dev>.175);
            plot(errev(1:end-1),diff(errev),'ro');
            
            % get the first mismatch
            errev1 = find(dev>.175,1,'first');
            errmag = dev(errev1);
                                
            if cnt > 4
                fprintf(['ERROR: makeSessionData.m line 332: correction ' ...
                         'count exceeded 4\n']);
                keyboard
            end
                 
            if ~isempty(errev1)

                if errev1 == length(bts) |  errev1 == length(bts)-1
                    % try fixing it by adjusting the behavioral event (this
                    % assumes neural events are ground truth
                    bts(errev1:end) = bts(errev1:end)+errmag;
                    
                elseif all(abs(bts((errev1+1:errev1+4)) - ...
                               nts((errev1:errev1+3))) < .1)
                    % if most of the remaining errors are low, this means
                    % that there is a missing neural event, so add
                    % one in
                    nts = [nts(1:errev1-1); bts(errev1); nts(errev1:end)];
                    errI(cnt) = errev1;
                    
                elseif diff(nts(errev1-1:errev1)) < 3
                    % for short IEI
                    bts(errev1:end) = bts(errev1:end)-errmag;
                    
                else
                    % try fixing it by adjusting the behavioral event (this
                    % assumes neural events are ground truth
                    bts(errev1:end) = bts(errev1:end)+errmag;
                    
                end

            end
            
            cnt = cnt + 1;
            drawnow;
            
        end
        
        % 3)
        % if events were added, append them to the stimOn vectors
        if ~isempty(errI)
            if ~any(round(stimOff - stimOn,3) ~= .005)
                % all events are .005
                stimOn = nts + stimOn(1);
                for i = 1:length(errI)
                    stimOff = [stimOff(1:errI(i)-1); stimOn(errI(i))+.005; ...
                               stimOff(errI(i):end)];
                end
            else
                keyboard
                % figure out if this is a target or onset event
                stimOn = nts + stimOn(1);
            end
        end
           
        
        
        %%
        % from here on depends on the event structure
        if ~any(round(stimOff - stimOn,3) ~= .005)
            %% if all of the events are .005, continue as usual:
            
            % some checks
            if length(stimOn)/2 ~= length(goodTrials)
                
                fprintf('Event mismatch: %5.2f stim onsets found, but %d trials.\n', ...
                        length(stimOn)/2,length(goodTrials));
                
                % check for one extra event from an aborted trial
                diffTrials = length(goodTrials) - length(stimOn)/2;
                if rem(diffTrials,1) ~= 0
                    
                    % if there is one extra event (stim onset probably from
                    % aborting a trial halfway through)
                    fprintf(['Luckily, one trial was aborted early in behavior, ' ...
                             'discarding that recording event!\n']);
                    stimOn = stimOn(1:end-1);
                    stimOff = stimOff(1:end-1);

                end
                
                %% fewer recording events
                if length(stimOn)/2 < length(goodTrials)
                    
                    % if there are fewer recording events than behavior events
                    fprintf(['It looks like the recording system missed some ' ...
                             'behavior events, let''s assume that the ' ...
                             'recording was stopped early, and remove those ' ...
                             'trials...\n']);
                    
                    % check for the number of trials missing
                    diffTrials = length(goodTrials) - length(stimOn)/2;
                    if rem(diffTrials,1) == 0
                        
                        % if a whole number, remove that number of trials
                        trialType(end-diffTrials+1:end,:) = [];
                        goodTrials(end-diffTrials+1:end) = [];
                        response(end-diffTrials+1:end) = [];
                        abort(end-diffTrials+1:end) = [];
                        
                    else
                        
                        keyboard
                        
                    end
                    
                end
                
                %% fewer behavior events
                if length(stimOn)/2 > length(goodTrials)
                    
                    diffTrials = length(goodTrials) - length(stimOn)/2;
                    fprintf(['Looks like there are extra recording computer ' ...
                             'events. Typically, this is impossible, but you ' ...
                             'should check to see what happened...\n']);
                    stimOn = stimOn(1:end-2);
                    stimOff = stimOff(1:end-2);
                    
                    if length(stimOn)/2 - length(goodTrials) > 10
                        keyboard
                    end
                    
                end
                
            elseif length(lickOn) ~= length(lickTimes)
                
                fprintf('Event mismatch: %d licks found, but %d licks recorded.\n', ...
                        length(lickOn),length(lickTimes));
                %sLicks = lickOn - lickOn(1);
                %bLicks = (lickTimes - lickTimes(1)) / 1e6;
                lickOn = lickOn(1:end-1);
                lickOff = lickOff(1:end-1);
                
            elseif length(rewdOn) ~= sum(response==1 & trialType(:,1)>0)
                
                fprintf('Event mismatch: %d rewards found, but %d rewards recorded.\n', ...
                        length(rewdOn),sum(response==1 & trialType(:,1)> ...
                                           0));
                
            end
            
            % stim and target indices
            stimI = [1:2:length(stimOn)];
            targI = [2:2:length(stimOn)];


        else
            %% if not all of the events are .005:
            
            % check for extra neural events
            if (length(stimOn) - length(bts)) == 1
                % remove last neural event
                stimOn(end) = [];
                stimOff(end) = [];
            end
            
            evd = round(stimOff - stimOn,3);
            stimI = find(evd == .005);
            targI = find(evd == .01);
            
            % check for event mismatch
            diffTrials = length(trialType) - length(stimI);
            if diffTrials > 0
                % check for extra behavior events (usually means recording ended early)
                fprintf(['Found %d behavioral trials, but %d neural ones...\n' ...
                         '(this usually means the recording was stopped ' ...
                         'before the behavior)\nRemoving %d extra behavior ' ...
                         'trials...\n'],length(trialType), ...
                        length(stimI),abs(diffTrials));
                trialType(end-diffTrials+1:end,:) = [];
                goodTrials(end-diffTrials+1:end) = [];
                response(end-diffTrials+1:end) = [];
                abort(end-diffTrials+1:end) = [];
                
            elseif diffTrials == -1
                % check for one extra neural trial
                fprintf(['Found %d behavioral trials, but %d neural ones, ' ...
                         'removing %d extra neural trial\n'],...
                        length(trialType),length(stimI),abs(diffTrials));
                
                % first check if the last trial had a target event too,
                % if so remove it, then remove the offending stimulus event
                dTargs = stimOn(targI) - stimOn(stimI(end));
                if any(dTargs > 0 & dTargs < 5)
                    targI(end) = [];
                end
                stimI(end) = [];
                
                
            elseif diffTrials < -1
                % error if more than one neural event missing
                fprintf(['%d extra neural events, this usually means ' ...
                         'something went wrong...\n'],abs(diffTrials));
                keyboard
                
                
            end
            
            
            % will probably need to do some checks here...
        end
        
        % remove the recording offset, unless there is a bug with the start time
        if mod(startTime/1e4,1) > 0
            stimOn = stimOn - startTime/fs;
            stimOff = stimOff - startTime/fs;
            lickOn = lickOn - startTime/fs;
            lickOff = lickOff - startTime/fs;
            rewdOn = rewdOn - startTime/fs;
            rewdOff = rewdOff - startTime/fs;
        end
        

        % get onset and target events
        noiseOn = stimOn(stimI);
        targOn = stimOn(targI);
        
        % find the index of trials which had targets
        for i = 1:length(noiseOn)
            dtarg = targOn - noiseOn(i);
            nAbort(i) = ~any(dtarg < 5 & dtarg > 0);
        end
        abortMismatch = find(nAbort' ~= abort);
        
        % if the recording stopped early, it may have missed the last
        % few trials, so there will be a mismatch between the behavior
        % and neural data. remove the last trials if this happens
        if ~isempty(abortMismatch)
            trialType(abortMismatch,:) = [];
            response(abortMismatch) = [];
            goodTrials(abortMismatch) = [];
            abort(abortMismatch) = [];
            noiseOn(abortMismatch) = [];
        end
        
        % if there is an extra stimulus onset relative to behavior,
        % remove the last neural trial:
        if (length(targOn) - length(trialType)) == 1
            % check if this trial had a target after
            targOn(end) = [];
        end
        if (length(noiseOn) - length(trialType)) == 1
            % check if this trial had a target after
            noiseOn(end) = [];
        end
        
        % assign target times to each trial (fill in nans where the
        % trial was aborted)
        targetOn = nan(size(noiseOn));
        targetOn(~abort) = targOn;
        
        
        
        %% process licks
        % tag lick events according to the trial they occurred in
        lickEdges = [noiseOn; inf];
        [~,~,index] = histcounts(lickOn,lickEdges);
        
        RT = [];
        lickTT = [];
        lickTrial = [];
        lickClockTime = [];
        respWindows = [];
        if ~isempty(index)
            % get the trial index and times of the first lick in each trial
            respWindows = [targetOn targetOn + b.params.respD];
            cnt = 0;
            for i = 1:length(respWindows)
                
                % get the first lick in this response window
                ind = find(lickOn > respWindows(i,1) & ...
                           lickOn < respWindows(i,2),1,'first');
                
                % if this trial has licks in the response window
                if ~isempty(ind)
                    
                    % save out variables for this trial
                    cnt = cnt + 1;
                    RT(cnt) = lickOn(ind) - respWindows(i,1);
                    lickTT(cnt) = trialType(i,1);
                    lickTrial(cnt) = i;
                    lickClockTime(cnt) = lickOn(ind);
                    
                end
                
            end
            
            % check for mismatches between neural and behavioral data
            if any(ismember(lickTrial',find(response>0)) == 0)
                
                % first check if there's an extra response trial in the
                % neural data, and remove it
                missedResponse = ismember(lickTrial',find(response>0));
                RT = RT(missedResponse);
                lickTT = lickTT(missedResponse);
                lickTrial = lickTrial(missedResponse);
                lickClockTime = lickClockTime(missedResponse);
                
            elseif any(ismember(find(response>0),lickTrial) == 0)
                
                % then check if there is an extra response trial in the
                % behavioral data
                a = find(response > 0);
                missedResponse = ismember(a,lickTrial);
                mr = a(~missedResponse);
                
                % remove from all behavioral data
                noiseOn(mr) = [];
                targetOn(mr) = [];
                trialType(mr,:) = [];
                response(mr,:) = [];
                goodTrials(mr,:) = [];
                abort(mr,:) = [];
                
            end
        end
        

        %%%% SAVE SESSION INFO %%%%
        %% parameters
        params.mouse = mouse;
        params.session = session;
        params.pathNeural = root;
        params.pathBehavior = broot;
        params.recordingStart = startTime;
        params.blockStart = blockStart;
        params.recordingMessages = msgtext;
        params.fs = fs;
        
        % plot title, etc.
        params.plot.titleString = [params.mouse '-' ...
                            params.session.session];
        params.plot.lw = 1.5;
        
        % plot colors for each contrast state target value
        n = length(b.params.targetDBShift);
        grade = linspace(.6,0,n);
        if b.params.sd(1) > b.params.sd(2)
            params.cond = 'hilo';
            params.plot.colors = [.3 .3 .3; grade' grade' ones(n,1)];
            params.plot.contrastColor = [1 0 0; 0 0 1];
        else
            params.cond = 'lohi';
            params.plot.colors = [.3 .3 .3; ones(n,1) grade' grade'];
            params.plot.contrastColor = [0 0 1; 1 0 0];
        end
        params.task = task;
        
        % condition numbers
        params.nLvl = length(unique(trialType(:,1)));
        params.nOff = length(unique(trialType(:,2)));
        
        % target offset (stimuli made before July 16, 2018 started 25ms
        % BEFORE the event)
        dt = datetime(session.session,'InputFormat','yyMMdd');
        if dt < datetime(2018,07,17)
            % targets come 25ms early
            params.offsets = b.params.noiseD - ...
                b.params.baseNoiseD - ...
                b.params.chordDuration;
            
        else
            % targets come on time
            params.offsets = b.params.noiseD - ...
                b.params.baseNoiseD;
            
        end
        
        % levels
        if size(b.params.targetDBShift,1)==2
            params.SNR = [[-inf; -inf] b.params.targetDBShift];
        else
            params.SNR = [-inf b.params.targetDBShift];
        end
        
        % other stim parameters
        params.stimInfo = b.params;
        
        % padding around trial offset and onset
        params.padding = .1;
        
        % session id
        params.sessionID = session.sessionInfo(end);
        
        %% events
        event.ts = ev.ts(ev.ts >= blockStart);
        event.state = ev.state(ev.ts >= blockStart);
        event.stim = [stimOn stimOff];
        event.lick = [lickOn lickOff];
        event.rewd = [rewdOn rewdOff];
        event.trialOn = noiseOn;
        event.targOn = targetOn;
        
        
        %% behavior
        behavior.task = task;
        behavior.matFile = session.behavior;
        behavior.txtFile = fn;
        behavior.trialInfo = t;
        behavior.eventTimes = ts;
        behavior.eventTags = tags;
        behavior.stimInfo = b.params;
        behavior.RT = RT;
        behavior.lickTT = lickTT;
        behavior.lickTrial = lickTrial;
        behavior.lickClockTime = lickClockTime;
        behavior.lickIndex = index;
        behavior.goodTrials = goodTrials;
        behavior.response = response;
        behavior.trialType = trialType;
        behavior.abort = abort;
        behavior.level = level;
        
        

        %% SPIKES
        % load spike info
        d = fullfile(root,'continuous','Rhythm_FPGA-100.0');
        recType = 'Rhythm_FPGA-100.0';
        if ~exist(d,'dir')
            % check for bandpassed data
            if exist(fullfile(root,'continuous',...
                              'Bandpass_Filter-107_100.0'),...
                     'dir');
                d = fullfile(root,'continuous','Bandpass_Filter-107_100.0');
                recType = 'Bandpass_Filter-107_100.0';
            else
                % for openEphys format, root is where spike data is saved
                d = root;
                recType = 'openephys';
            end

        end
        SPIKES = double(readNPY([d filesep 'spike_times.npy'])) / fs;
        clusters = double(readNPY([d filesep 'spike_clusters.npy']));
        
        % check for cluster groups
        if ~exist([d filesep 'cluster_groups.csv'])
            % for .tsv and .csv files
            if exist([d filesep 'cluster_group.tsv'])
                fn = [d filesep 'cluster_group.tsv'];
                
            else
                % if no manually sorted files, check for the raw
                % kilosort output
                fprintf('No sorted data was found!\n');
                if exist([d filesep 'cluster_KSLabel.tsv'])
                    % option to use that data
                    flag = [];
                    while isempty(flag)
                        flag = input(['Found some raw KiloSorted data; press ' ...
                                      '[y] to use, or [n] to cancel: '],'s');
                    end
                    
                    if strcmp(flag,'y') | strcmp(flag,'Y')
                        fprintf(['Continuing with raw KiloSort ' ...
                                 'output...\n']);
                        fn = [d filesep 'cluster_KSLabel.tsv'];
                        
                    else
                        fprintf('No sorted data was found!\n');
                        keyboard
                        
                    end
                end
            end
            % load tab delimited data
            tmp = importdata(fn,'\t');

        else
            % load comma delimited data
            fn = [d filesep 'cluster_groups.csv'];
            tmp = importdata(fn,',');
            
        end
        
        % extract labels
        labels = vertcat(tmp(2:end));
        if isempty(labels)
            fprintf('No labels found in %s... stopping.\n',fn);
            keyboard
        end
        
        %%%% ORGANIZE SPIKES %%%%
        cnt = 0;
        for c = 1:length(labels)
            
            % get cluster ID and quality
            str = strsplit(labels{c});
            ID(c) = str2num(str{1});
            qstr{c} = str{2};
            switch qstr{c}
              case 'noise'
                q(c) = 0;
              case 'good'
                q(c) = 1;
                lbl = 'su';
              case 'mua'
                q(c) = 2;
                lbl = 'mu';
            end
            
            % if the cluster is noise, skip it
            fprintf('\tCluster %02d/%02d...',c,length(labels));
            if q(c) == 0
                fprintf(' skipping...\n');
                continue;
            end
            cnt = cnt + 1;
            
            % extract spikes for this cell
            spks = SPIKES(clusters == ID(c));
            meanFR = length(spks)/(spks(end)-spks(1));
            fprintf('FR = %04.2fHz\n',meanFR);
                        
            % write out cell info
            cellID = sprintf('%s_%d_%02d_%03d_%s',mouse,session.sessionInfo(end),...
                             cnt,ID(c),lbl);
            cellInfo(cnt,:) = {mouse, session.sessionInfo(end), cnt, ID(c), lbl, ...
                               q(c), cellID, meanFR, ...
                               session.behavior, session.neural, session.type};
            
            % save out all spikes to a cell array
            spikes{cnt} = spks;
        end
        params.nCells = cnt;
        
        
        
        %% ADD WAVEFORMS
        waveFile = ...
            dir(fullfile(params.pathNeural,'continuous', ...
                         recType, ...
                         'mean_waveforms.mat'));
        
        if ~isempty(waveFile)
            
            fprintf('\tAdding waveforms!\n'); 
            
            % load waveforms
            wv = load(fullfile(waveFile.folder,waveFile.name));
                        
            % unit index
            I = q>0;
            
            % add indexed vars to data structure
            if isfield(wv,'chanMap')
                data.wave.chanMap = wv.chanMap;
            else
                data.wave.chanMap =[];
            end
            data.wave.waveform = wv.mw(I,:,:);
            data.wave.max = wv.mx(I,:);
            data.wave.noise_level = wv.noise_level;
            data.wave.percentile95 = wv.s95(I,:,:,:);
            data.wave.spkcount = wv.sn(I)';
            data.wave.snr = wv.snr(I,:);
            data.wave.stddev = wv.sw(I,:,:);
            data.wave.peak = max(abs(data.wave.waveform),[],3);
            [~,data.wave.peakchan] = max(data.wave.peak,[],2);
            for i = 1:size(data.wave.waveform,1)
                
                % waveform
                mxwave = squeeze(data.wave.waveform(i,data.wave.peakchan(i),:));
                
                % interpolate
                t0 = (1:length(mxwave)) / params.fs;
                t1 = linspace(t0(1),t0(end),1000);
                mxwave1 = interp1(t0,mxwave,t1,'cubic');
                
                % get the trough, post-spike peak, post-spike inflection
                [trVal,trI] = min(mxwave1);
                tmp = mxwave1; tmp(1:trI) = nan;
                [pkVal,pkI] = max(tmp);
                tmp = mxwave1; tmp(1:pkI) = nan;
                [infVal,infI] = min(abs(tmp-0));
                
                % full width half min
                hm = trVal / 2;
                tmp = mxwave1-hm; tmp(1:trI) = nan;
                [~,hI(1)] = min(abs(tmp-0));
                tmp = mxwave1-hm; tmp(trI:end) = nan;
                [~,hI(2)] = min(abs(tmp-0));
                
                data.wave.peakwaveform(i,:) = mxwave;
                data.wave.trough_peak(i,1) = (t1(pkI) - t1(trI));
                data.wave.peak_inflect(i,1) = (t1(infI) - t1(pkI));
                data.wave.FWHM(i,1) = abs(diff(t1(hI(1:2))));
                data.wave.fs = params.fs;
                
            end
            
            
            
            % plot results
            subplot(1,5,4); 
            [~,pksort] = sort(data.wave.peakchan);
            imagesc(data.wave.peak(pksort,:))
            
            subplot(1,5,5); plot(data.wave.peakwaveform');
                        
            
        else
            
            nu = sum(q>0);
                        
            fprintf('\tNo waveforms...\n');
            data.wave.chanMap = [];
            data.wave.waveform = [];
            data.wave.max = [];
            data.wave.noise_level = [];
            data.wave.percentile95 = [];
            data.wave.spkcount = nan(nu,1);
            data.wave.snr = [];
            data.wave.stddev = [];
            data.wave.peak = [];
            data.wave.peakchan = nan(nu,1);
            data.wave.peakwaveform = nan(nu,63);
            data.wave.trough_peak = nan(nu,1);
            data.wave.peak_inflect = nan(nu,1);
            data.wave.FWHM = nan(nu,1);
            data.wave.fs = nan(1,1);

            
        end
        
        % save the alignment figure
        fn = sprintf('~/data/kilosort/_event_alignment_plots/%s-%s.pdf',...
                     mouse,session.session);
        saveFigPDF(f123,[1200 300],fn);
        
        % save out data
        fprintf(sprintf('\tSaving to %s...\n',output));
        data.cellInfo = cellInfo;
        data.events = event;
        data.behavior = behavior;
        data.params = params;
        data.params.fnout = fnout;
        sessionInfo = params;
        save(fnout, 'data','sessionInfo','spikes');
        
    else
        if nodataflag
            fprintf('WARNING: No MUAs or SUs found... skipping\n');
        else
            fprintf('WARNING: NO CURATED DATA FOUND... skipping\n');
        end

        fnout = [];
        
    end

else
    fprintf('File %s already exists...\n',fnout);
end




















%% DIGRESSION: Matching recording and behavior system times
%
% *** Basic conclusion: the arduino clock drifts a lot! we
% shouldn't trust the timestamps we get very much...

% attempt to match up neural and behavioral times
% plots below indicate quite a bit of drift between the two clocks
% that is not constant :/
%plot((lickOn-stimOn(1))-(lickTimes/1e6))
%plot(lickOn-stimOn(1),lickTimes / 1e6)
%plot(lickOn-stimOn(1),ones(size(lickOn)),'ko')
%plot(lickTimes / 1e6,ones(size(lickTimes)) * 2,'bo')

% try matching them per trial
% running this indicates some pretty severe clock drift on the
% arduino (~3-5ms over a 4-7s trial), but as long as everything
% matches up with the recording system, we can trust that clock for
% aligning to our spike data
%  for i = 1:length(noiseOn)-1
%      % for each stimulus onset, find licks, rewards and check timing
%      lickEvs{i} = lickOn(lickOn > noiseOn(i) & lickOn < noiseOn(i+1));
%      lickEvs{i} = lickEvs{i} - noiseOn(i);
%      rewdEvs{i} = rewdOn(rewdOn > noiseOn(i) & rewdOn < noiseOn(i+1));
%      rewdEvs{i} = rewdEvs{i} - noiseOn(i);
%      
%      dt = (t(i).lickTimes - t(i).stimTimes)/1e6 - lickEvs{i};
%      disp(dt)
%      pause
%  end
%
%% ADDENDUM TO DIGRESSION
% turns out that the drift was due to sloppy plotting in
% matlab, vectorizing the plots made them much quicker, and so
% there is no longer any added delay as the number of trials to
% plot increases
