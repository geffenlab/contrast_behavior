function [res,r] = run_choice(spikeData,sessionData,ops)

% function res = run_choice(spks,cellinfo,sessionData,ops)

sig_cells = ops.sig_neurons;
resFile = fullfile(ops.resDir,ops.resFile);

% get neurons in the task that have good waveforms
included_cells = ops.include & contains(spikeData.cellinfo(:,end),ops.task);
spikes = spikeData.spikes(included_cells);
cellInfo = spikeData.cellinfo(included_cells,:);


for c = 1:length(spikes)
    
    fprintf('Neuron %d/%d: ',c,length(spikes)); tic;
    fn = fullfile(ops.resDir,'_choice',[cellInfo{c,7} '.mat']);
    
    if ~exist(fn)
        fprintf('computing cp... ');
        
        % find session info for this neuron
        sI = find(cat(1,sessionData.session.sessionID) == cellInfo{c,2});
        session = sessionData.session(sI);
        events = sessionData.events(sI);
        behavior = sessionData.behavior(sI);
        
        % compute PSTH relative to target onset
        [PSTH,raster,trials] = makePSTH(spikes{c},events.targOn,ops.edges);
        
        
        % index for hits and misses during target trials
        hi = behavior.trialType(:,1) > 0 & behavior.response == 1;
        mi = behavior.trialType(:,1) > 0 & behavior.response == 0;
        
        % clear
        clear nhits nmiss cp cp_pct cp_sig;
        
        % choice probability
        uVols = unique(behavior.trialType(behavior.trialType(:,1)>0,1));
        for i = 1:length(uVols)
            
            hitI = behavior.trialType(:,1) == uVols(i) & hi;
            missI = behavior.trialType(:,1) == uVols(i) & mi;
            nhits(i) = sum(hitI);
            nmiss(i) = sum(missI);
            
            % compute cp over time for each volume
            for t = 1:length(ops.timeCent)
                
                % index
                timeI = ops.time > ops.timeCent(t)-ops.timeWindow & ...
                        ops.time <= ops.timeCent(t)+ops.timeWindow;

                
                if sum(hitI) > 1 & sum(missI) > 1
                    % average spike rate
                    hit_spks = mean(PSTH(hitI,timeI),2);
                    miss_spks = mean(PSTH(missI,timeI),2);
                    
                    [cp(i,t),~,~,~,cp_pct(i,t,:),~,~,cp_sig(i,t)] = ...
                        bootROC(miss_spks,hit_spks,[],100);
                else
                    cp(i,t) = nan;
                    cp_pct(i,t,:) = [nan nan];
                    cp_sig(i,t) = false;
                end
            end
            
            
        end
        
        % results
        r.cellID = cellInfo{c,7};
        r.cellInfo = cellInfo(c,:);
        r.sessionID = cellInfo{c,2};
        r.contrast = session.cond;
        r.cp = cp;
        r.cp_pct = cp_pct;
        r.cp_sig = cp_sig;
        save(fn,'r','ops');
        
    else
        fprintf('file found, loading... ');
        load(fn,'r');
        
    end
    
    toc;
    
    res(c) = r;
    
end





