function [res,r,ops] = run_offset(spikeData,sessionData,ops)

sig_cells = ops.sig_neurons;
resFile = ops.resFile;

% get neurons in the psychometric task that have good waveforms
included_cells = ops.include & contains(spikeData.cellinfo(:,end),'offset');
spikes = spikeData.spikes(included_cells);
cellInfo = spikeData.cellinfo(included_cells,:);

figure;

%% SINGLE NEURONS
if ~exist(fullfile(ops.resDir,resFile),'file')
    

    % first run analysis of single neurons
    t0 = tic;
    fprintf('OFFSET ANALYSIS: single neurons\n');
    for i = 1:length(spikes)
        
        fprintf('\tCell %d/%d ',i,length(spikes));
        
        % find session info, events, behavior
        s = sessionData.session(cellInfo{i,2} == [sessionData.session.sessionID]);
        e = sessionData.events(cellInfo{i,2} == [sessionData.session.sessionID]);
        b = sessionData.behavior(cellInfo{i,2} == [sessionData.session.sessionID]);
        
        % get PSTH
        [PSTH,raster,trials,PSTHs] = makePSTH(spikes{i},e.targOn,ops.edges,ops.smooth);
        t_spikes = mean(PSTH(:,ops.time > ops.target(1) & ops.time <= ops.target(2)),2); % not including <= (equals) can induce rounding errors???
        [~,sortI] = sortrows(b.trialType);
        
        % trials to include
        include = ~b.abort & b.goodTrials;
        
        % unique trials
        uVols = 0:s.nLvl-1;
        uOffs = 1:s.nOff;
        cnt = 1; uTrials = [];
        for j = 1:length(uVols)
            for k = 1:length(uOffs)
                uTrials(cnt,1) = uVols(j);
                uTrials(cnt,2) = uOffs(k);
                cnt = cnt + 1;
            end
        end
                
        % mean rate per condition
        t_spikes_snr = nan(s.nLvl,s.nOff);
        auc = nan(s.nLvl-1,s.nOff);
        auc_pct = nan(s.nLvl-1,s.nOff,2);
        auc_sig = nan(s.nLvl-1,s.nOff);
        for j = 1:s.nLvl
            for k = 1:s.nOff
                
                % trial index
                I = b.trialType(:,1) == (j-1) & b.trialType(:,2) == (k) & ...
                    include;
                
                % mean spike rate and behavior detection rate
                t_spikes_snr(j,k) = mean(t_spikes(I));
                
                % auc
                if j == 1
                    ndist = t_spikes(I);
                else
                    sdist = t_spikes(I);
                    [auc(j-1,k),~,~,~,auc_pct(j-1,k,:),~,~,auc_sig(j-1,k)] = ...
                        bootROC(ndist,sdist,[],500);
                end
            end
        end
        
        % compute behavioral performance
        tmp = compute_yesno_perf(b.response(include),b.trialType(include,[1 2]),uTrials,...
                                 true);
        br_adj = reshape(tmp,s.nOff,s.nLvl)';
        tmp = compute_yesno_perf(b.response(include),b.trialType(include,[1 2]),uTrials,...
                                 false);
        br = reshape(tmp,s.nOff,s.nLvl)';
        
        
        % single cell criterion classifier
        [pred, crit, rule] = crit_classifier(t_spikes(include),b.trialType(include,[1 2]));
        tmp = compute_yesno_perf(pred,b.trialType(include,[1 2]),uTrials,...
                                 true);
        critp_adj = reshape(tmp,s.nOff,s.nLvl)';
        tmp = compute_yesno_perf(pred,b.trialType(include,[1 2]),uTrials,...
                                 false);
        critp = reshape(tmp,s.nOff,s.nLvl)';
        
        % runtime plot
        t(i) = toc(t0); toc(t0);
        plot(diff(t));
        xlabel('Cell'); ylabel('Time (s)');
        xlim([0 length(spikes)]);
        title('run_offset.m progress');
        set(gca, 'yscale', 'log')
        drawnow;
        
        % save results
        res.single_cell(i).cellID = cellInfo{i,7};
        res.single_cell(i).PSTH = PSTH;
        res.single_cell(i).PSTHs = PSTHs;
        res.single_cell(i).trial_fr = t_spikes;
        res.single_cell(i).brate = br;
        res.single_cell(i).brate_adj = br_adj;
        res.single_cell(i).fr = t_spikes_snr;
        res.single_cell(i).auc = auc;
        res.single_cell(i).auc_pct = auc_pct;
        res.single_cell(i).sig = auc_sig;
        res.single_cell(i).critp = critp;
        res.single_cell(i).critp_adj = critp_adj;
        res.single_cell(i).vols = s.SNR;
        res.single_cell(i).offs = s.offsets;
        
        % f = plot_single_cell_offset(cellInfo{i,7},res.single_cell(i),spikeData,sessionData,ops);
        
        
        
    end
    
    % save
    save(fullfile(ops.resDir,resFile),'res','ops','-v7.3');
    
else
    
    fprintf('Found res file... loading ...\n');
    
    % load
    load(fullfile(ops.resDir,resFile));
    
end

if ~isfield(res,'pop')

    % consider only psych sessions
    s = sessionData.session(contains({sessionData.session.task},'offset'));
    b = sessionData.behavior(contains({sessionData.session.task},'offset'));
    e = sessionData.events(contains({sessionData.session.task},'offset'));

    % number of significant neurons per recording
    sig = {res.single_cell.sig};
    sig = cellfun(@(x) sum(x,2),sig,'UniformOutput',false);
    nsig = [sig{:}];

    % significant cells defined as having at least one significant
    % response at low and at high volume
    sigCells = all(nsig>0)';

    i = 0;
    for sess = 1:length(s)
        
        % cells to include       
        cellI = ([cellInfo{:,2}] == s(sess).sessionID)';
        if sig_cells
            fprintf('POPULATION RESULTS USING ONLY SIG NEURONS\n');
            cellI = cellI & sigCells;
            plot_sig = 'sig_';
        end
        ops.sig_neurons = sig_cells;
        
        % if there are 3 or more cells in this session
        if sum(cellI) >= 3
            
            
            %% analysis
            fprintf('SVM session %d/%d... ',sess,length(s)); tic;
            
            i = i + 1;
            
            % include good trials
            include = ~b(sess).abort & b(sess).goodTrials;
            
            % behavioral performance
            br = mean(cat(3,res.single_cell(cellI).brate),3);
            br_adj = mean(cat(3,res.single_cell(cellI).brate_adj),3);
            
            % volumes
            vols = b(sess).trialType(include,1);
            
            % offsets
            offs = b(sess).trialType(include,2);
            
            % unique trials
            uVols = 0:s(sess).nLvl-1;
            uOffs = 1:s(sess).nOff;
            cnt = 1; uTrials = [];
            for j = 1:length(uVols)
                for k = 1:length(uOffs)
                    uTrials(cnt,1) = uVols(j);
                    uTrials(cnt,2) = uOffs(k);
                    cnt = cnt + 1;
                end
            end
            
            % target response matrix
            X = [res.single_cell(cellI).trial_fr];
            X = X(include,:);
            
            % labels (target present or not)
            labels = vols > 0;
            
            % dummy code unique trial combos
            grps = num2str([vols offs]);
            
            % run svm
            cvp = cvpartition(grps,'KFold',10,'Stratify',true);
            mdl = fitcsvm(X,labels,'KernelFunction','linear','CVPartition',cvp,...
                          'Standardize',true,...
                          'Prior',[sum(labels == 0) / length(labels) ...
                                sum(labels == 1) / length(labels)]);
            elabel = kfoldPredict(mdl);
            
            % train CD
            y = vols;
            [projection cd wt] = trainCD(X,y,[],ops);
            normp = normalize(projection,'range');
            
            % criterion classifier using the projection
            [pred, crit, rule] = crit_classifier(normp,vols);
            tmp = compute_yesno_perf(pred,[vols offs],uTrials,...
                                           true);
            critr_adj = reshape(tmp,s(sess).nOff,s(sess).nLvl)';
            tmp = compute_yesno_perf(pred,[vols offs],uTrials,...
                                       false);
            critr = reshape(tmp,s(sess).nOff,s(sess).nLvl)';
            
            
            % svm performance
            tmp = compute_yesno_perf(elabel,[vols offs],uTrials,...
                                           true);
            psvm_adj = reshape(tmp,s(sess).nOff,s(sess).nLvl)';
            tmp = compute_yesno_perf(elabel,[vols offs],uTrials,...
                                          false);
            psvm = reshape(tmp,s(sess).nOff,s(sess).nLvl)';

            
            % population performance per snr (svm and AUC)
            pauc = nan(length(uVols)-1,length(uOffs));
            pauc_sig = nan(length(uVols)-1,length(uOffs));
            pauc_pct = nan(length(uVols)-1,length(uOffs),2);
            for j = 1:length(uVols)
                for k = 1:length(uOffs)
                    I = vols == uVols(j) & offs == uOffs(k);
                                        
                    % auc
                    if j == 1
                        ndist = normp(I);
                    else
                        sdist = normp(I);
                        [pauc(j-1,k),~,~,~,pauc_pct(j-1,k,:),~,~,pauc_sig(j-1,k)] = ...
                            bootROC(ndist,sdist,[],500);
                    end
                end
            end
            
            % outputs
            res.pop(i).beh_rate = br;
            res.pop(i).beh_rate_adj = br_adj;
            res.pop(i).svm_rate = psvm;
            res.pop(i).svm_rate_adj = psvm_adj;
            res.pop(i).critp = critr;
            res.pop(i).critp_adj = critr_adj;
            res.pop(i).auc = pauc;
            res.pop(i).auc_pc = pauc_pct;
            res.pop(i).auc_sig = pauc_sig;
            res.pop(i).projection = normp;
            res.pop(i).mean_auc = mean(...
                cat(3,res.single_cell(cellI&sigCells).auc),3);
            res.pop(i).mean_fr = mean(...
                cat(3,res.single_cell(cellI&sigCells).fr),3);
            res.pop(i).mean_critp = mean(...
                cat(3,res.single_cell(cellI&sigCells).critp),3);
            res.pop(i).mean_critp_adj = mean(...
                cat(3,res.single_cell(cellI&sigCells).critp_adj),3);
            res.pop(i).snr = s(sess).SNR;
            res.pop(i).offs = s(sess).offsets;
            res.pop(i).sessID = s(sess).sessionID;
            
            toc
            
            
            
            %% plotting
            plotOn = true;
            if plotOn
                
                f1 = figure(1122); set(gcf,'Position',[0 0 650 500]); clf;
                plot_session_summary_offset(s(sess).sessionID,res,spikeData,sessionData,ops);
                
                fn = sprintf('./_plots/_session/offset_%s_%d.pdf',s(sess).mouse,s(sess).sessionID);
                saveFigPDF(f1,[650 500],fn);
                
            end
            
        end
        
        
    end

    save(fullfile(ops.resDir,resFile),'res','ops','-v7.3');
    
else
    
    load(fullfile(ops.resDir,resFile));
    
end



%% summary stats
if ~exist('r','var')
    
    fprintf('Post-processing...\n');

    % restructure to match different offset times, first get unique
    % offsets
    uOffs = unique(round([res.pop(:).offs],3));

    % set fields to use
    fields = {'beh_rate','svm_rate','critp','auc','auc_sig'...
              'mean_auc','mean_fr','mean_critp',...
              'beh_rate_adj','svm_rate_adj','critp_adj','mean_critp_adj'};

    % build nan padded data
    clear r;
    for i = 1:length(res.pop)
        
        this_offset = round(res.pop(i).offs,3);
        offI = ismember(uOffs,this_offset);
        
        % add offsets
        r.offs(i,:) = uOffs;
        
        % count number of significant neurons/total neurons (first need
        % to index into the single_cell res structure
        sI = contains({res.single_cell.cellID}, ...
                      num2str(res.pop(i).sessID));
        sig_neurons = sum(cat(3,res.single_cell(sI).sig),3,'omitnan');
        r.tot_neurons(i,1) = sum(sI);
        r.sig_neurons(i,:,:) = nan(2,length(uOffs));
        r.sig_neurons(i,1,offI) = sig_neurons(1,:);
        r.sig_neurons(i,2,offI) = sig_neurons(2,:);
        r.prop_sig_neurons(i,:,:) = r.sig_neurons(i,:,:)./r.tot_neurons(i);
        
        
        % add target fields
        for j = 1:length(fields)
            sz = size(res.pop(i).(fields{j}),1);
            if sz == 0
                sz = size(res.pop(1).(fields{j}),1);
            end
            r.(fields{j})(i,:,:) = nan(sz,length(uOffs));
            for k = 1:size(res.pop(i).(fields{j}),1)
                r.(fields{j})(i,k,offI) =  res.pop(i).(fields{j})(k,:);
            end
        end
    end

    % add in other useful info
    r.sessionID = vertcat(res.pop.sessID);
    for i = 1:length(res.pop)
        
        % mouse and contrast
        sI = ismember([sessionData.session.sessionID],r.sessionID(i));
        r.mouse{i,1} = sessionData.session(sI).mouse;
        r.contrastI(i,1) = strcmp(sessionData.session(sI).cond,'lohi');
        
    end
    r.fa = mean(r.(fields{1})(:,1,:),[2 3],'omitnan');

    % save formatted data and new options
    save(fullfile(ops.resDir,resFile),'res','r','ops','-v7.3');

else
    load(fullfile(ops.resDir,resFile));
    
end








































if false
    % plots for each mouse
    uM = unique(r.mouse,'row');
    for i = 1:length(uM)
        
        figure(i); clf;
        
        for j = 1:2
            
            I = contains(r.mouse,uM{i}) & r.contrastI == (j-1);
            
            if sum(I) > 0
                
                for k = 1:length(fields)
                    tmp = r.(fields{k})(I,:,:);
                    sz = size(tmp);
                    
                    subplot(length(fields),2,2*(k-1)+j); hold on;
                    
                    for ii = 1:sz(2)
                        % plot each line
                        if sz(2) < 3
                            colI = ii + 1;
                        else 
                            colI = ii;
                        end
                        
                        %reshape(tmp,[sz([1 3]) 2])
                        
                        plot(r.offs(I,:),squeeze(tmp(:,ii,:)),...
                             '.','color',colors{j}{colI});
                        x = mean(r.offs(I,:),1,'omitnan');
                        y = mean(squeeze(tmp(:,ii,:)),1,'omitnan');
                        plot(x(~isnan(y)),y(~isnan(y)),'color',colors{j}{colI});
                        
                    end
                    set(gca,'xtick',x(~isnan(y))); set(gca,'xticklabels', ...
                                                           num2str(x(~isnan(y))'* ...
                                                                   1000));
                    xlim([-.05 1.05]);
                    if k == length(fields)
                        xlabel('Time (s from switch)');
                    end
                    plot([0 0],ylim,'k');
                    ylabel(sprintf('%s',fields{k}),'interpreter', ...
                           'none');
                    plotPrefs;
                    if k == 1
                        title(uM{i});
                    end
                    
                end
            end
        end
    end
end

                
                
            
            
    
    
    

    



