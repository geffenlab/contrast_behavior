function [res,r] = run_psych(spikeData,sessionData,ops)

% function res = run_psych(spks,cellinfo,sessionData,ops)
sig_cells = ops.sig_neurons;
resFile = ops.resFile;

% get neurons in the psychometric task that have good waveforms
included_cells = ops.include & contains(spikeData.cellinfo(:,end),'psychometric');
spikes = spikeData.spikes(included_cells);
cellInfo = spikeData.cellinfo(included_cells,:);

%a = load('/Users/chris/chris-lab/projects/contrast_behavior/_data_back/res_psych_sigcells.mat');


%% SINGLE NEURONS
if ~exist(fullfile(ops.resDir,resFile),'file')
    
    % set up figure for internal use
    fh = figure(1234); clf;
    set(fh,'Visible','on');  

    % first run analysis of single neurons
    t0 = tic;
    fprintf('PSYCH ANALYSIS: single neurons\n');
    for i = 1:length(spikes)
        
        fprintf('\tCell %d/%d ',i,length(spikes));
        
        % find session info, events, behavior
        s = sessionData.session(cellInfo{i,2} == [sessionData.session.sessionID]);
        e = sessionData.events(cellInfo{i,2} == [sessionData.session.sessionID]);
        b = sessionData.behavior(cellInfo{i,2} == [sessionData.session.sessionID]);
        
        % get PSTH
        [PSTH,raster,trials,PSTHs] = makePSTH(spikes{i},e.targOn,ops.edges,ops.smooth);
        t_spikes = mean(PSTH(:,ops.time > ops.target(1) & ops.time <= ops.target(2)),2); % not including <= (equals) can induce rounding errors???     
        
        % trials to include
        include = ~b.abort & b.goodTrials;
        
        % mean rate per condition
        t_spikes_snr = nan(size(s.SNR));
        b_presp_snr = nan(size(s.SNR));
        auc = nan(1,size(s.SNR,2)-1);
        auc_sig = nan(1,size(s.SNR,2)-1);
        for j = 1:length(s.SNR)
            
            % trial index
            I = b.trialType(:,1) == (j-1) & include;
            
            % mean spike rate and behavior detection rate
            t_spikes_snr(j) = mean(t_spikes(I));
            
            % auc
            if j == 1
                ndist = t_spikes(I);
            else
                sdist = t_spikes(I);
                [auc(j-1),~,~,~,auc_pct(j-1,:),~,~,asig] = bootROC(ndist,sdist,[],500);
                auc_sig(j-1) = ~(.5 >= auc_pct(j-1,1) & .5 <= auc_pct(j-1,2));
            
            end
            
        end
        
        % compute behavioral performance
        uTrials = [0:s.nLvl-1]';
        br_adj = compute_yesno_perf(...
            b.response(include),b.trialType(include,1),uTrials,...
            true); 
        br = compute_yesno_perf(...
            b.response(include),b.trialType(include,1),uTrials,...
            false);   
        
        % single cell criterion classifier
        [pred, crit, rule] = crit_classifier(t_spikes(include),b.trialType(include,1));
        critr_adj = compute_yesno_perf(pred,b.trialType(include,1),uTrials,...
                                       true);
        critr = compute_yesno_perf(pred,b.trialType(include,1),uTrials,...
                                   false);
        
        % runtime plot
        t(i) = toc(t0); toc(t0);
        plot(diff(t));
        xlabel('Cell'); ylabel('Time (s)');
        xlim([0 length(spikes)]);
        set(gca, 'yscale', 'log')
        title('run_psych.m progress');
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
        res.single_cell(i).critp = critr;
        res.single_cell(i).critp_adj = critr_adj;
        res.single_cell(i).vols = s.SNR;
        
        % f = plot_single_cell_psych(cellInfo{i,7},res,spikeData,sessionData,ops);
        
    end
    
    % save res file
    save(fullfile(ops.resDir,resFile),'res','ops','-v7.3');

else
    
    fprintf('Res file found, loading... \n');
    
    % load res file
    load(fullfile(ops.resDir,resFile));
    
end

% use frozen auc percentiles for stable significance over multiple runs
tmp = load('./_data/psych_sig_cells.mat');
for i = 1:length(res.single_cell)
    res.single_cell(i).auc_pct = tmp.auc_pct{i};
    res.single_cell(i).sig = tmp.auc_sig{i};
end



if ~isfield(res,'pop')
    
    % set up figure for internal use
    fh = figure(1234); clf;
    set(fh,'Visible','on');
    
    %% population analysis per session
    % consider only psych sessions
    s = sessionData.session(contains({sessionData.session.task},'psychometric'));
    b = sessionData.behavior(contains({sessionData.session.task},'psychometric'));
    e = sessionData.events(contains({sessionData.session.task},'psychometric'));

    % significant cells defined as cells with significant responses
    % to two or more volumes
    sig = vertcat(res.single_cell.sig);
    sigCells = sum(sig,2) >= 2;

    i = 0;
    for sess = 1:length(s)
        
        plot_sig = [];
        
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
            br = mean(vertcat(res.single_cell(cellI).brate),1);
            br_adj = mean(vertcat(res.single_cell(cellI).brate_adj),1);
            
            % volumes
            vols = b(sess).trialType(include,1);
            
            % target response matrix
            X = [res.single_cell(cellI).trial_fr];
            X = X(include,:);     
            
            % labels
            labels = vols > 0;
            
            % run svm
            cvp = cvpartition(vols,'KFold',10,'Stratify',true);
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
            uTrials = [0:s(sess).nLvl-1]';
            [pred, crit, rule] = crit_classifier(normp,vols);
            critr_adj = compute_yesno_perf(pred,vols,uTrials,...
                                           true);
            critr = compute_yesno_perf(pred,vols,uTrials,...
                                       false);
            
            % svm performance
            psvm_adj = compute_yesno_perf(elabel,vols,uTrials,...
                                           true);
            psvm = compute_yesno_perf(elabel,vols,uTrials,...
                                          false);
            
            % performance per snr
            pauc = nan(1,length(uTrials)-1);
            pauc_pct = nan(length(uTrials)-1,2);
            pauc_sig = nan(1,length(uTrials));
            for j = 1:length(uTrials)
                I = vols == uTrials(j);
                
                if sum(I) > 0
                    % auc
                    if j == 1
                        ndist = normp(I);
                    else
                        sdist = normp(I);
                        [pauc(j-1),~,~,~,pauc_pct(j-1,:),~,~,pauc_sig(j-1)] = ...
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
                cat(1,res.single_cell(cellI&sigCells).auc),1);
            res.pop(i).mean_fr = mean(...
                cat(1,res.single_cell(cellI&sigCells).fr),1);
            res.pop(i).mean_critp = mean(...
                cat(1,res.single_cell(cellI&sigCells).critp),1);
            res.pop(i).mean_critp_adj = mean(...
                cat(1,res.single_cell(cellI&sigCells).critp_adj),1);
            res.pop(i).snr = s(sess).SNR;
            res.pop(i).sessID = s(sess).sessionID;
                        
            toc
            
                        
            %% plotting
            plotOn = true;
            if plotOn
                
                clf;
                plot_session_summary_psych(s(sess).sessionID,res,spikeData,sessionData,ops);

                if ~exist('./_plots/_session')
                    mkdir('./_plots/_session');
                end
                
                fn = sprintf('./_plots/_session/psych_%s%s_%d.pdf',...
                             plot_sig,s(sess).mouse,s(sess).sessionID);
                saveFigPDF(fh,[650 500],fn);
                
            end
            
        end
        
    end
    
    save(fullfile(ops.resDir,resFile),'res','ops','-v7.3');
    
else
    
    load(fullfile(ops.resDir,resFile));
    
end

% use frozen population auc percentiles for stable significance
% values over multiple runs
tmp = load('./_data/psych_sig_cells.mat');
for i = 1:length(res.pop)
    res.pop(i).auc_pc = tmp.pop_auc_pct{i};
    res.pop(i).auc_sig = tmp.pop_auc_sig{i};
end

if ~exist('r','var');
    
    fprintf('Post-processing (this can take a while)...\n');

    % build nanpadded data for unique volumes
    uVol = unique(cat(1,res.pop.snr));
    uVol_nn = uVol(2:end);

    % set fields to use
    fields = {'beh_rate','svm_rate','critp','auc','auc_sig'...
              'mean_auc','mean_fr','mean_critp',...
              'beh_rate_adj','svm_rate_adj','critp_adj','mean_critp_adj'};


    clear r;
    for i = 1:length(res.pop)
        
        fprintf('\tSession %d/%d... ',i,length(res.pop)); tic;
        
        this_vol = round(res.pop(i).snr,3);
        volI = ismember(uVol,this_vol);
        volI_nn = ismember(uVol_nn,this_vol);
        
        % add offsets
        r.vols(i,:) = uVol;
        r.vols_nn(i,:) = uVol_nn;
        
        % count number of significant neurons/total neurons (first need
        % to index into the single_cell res structure
        sI = contains({res.single_cell.cellID}, ...
                      num2str(res.pop(i).sessID));
        sig_neurons = sum(cat(1,res.single_cell(sI).sig),'omitnan');
        r.tot_neurons(i,1) = sum(sI);
        r.sig_neurons(i,:) = nan(1,length(uVol_nn));
        r.sig_neurons(i,volI_nn) = sig_neurons;
        r.prop_sig_neurons(i,:) = r.sig_neurons(i,:)./r.tot_neurons(i);
        
        
        % add target fields
        for j = 1:length(fields)
            sz1 = size(res.pop(i).(fields{j}),1);
            sz2 = size(res.pop(i).(fields{j}),2);
            
            if sz2 == 7
                vI = volI;
            elseif sz2 == 6
                vI = volI_nn;
            elseif sz2 == 0
                if length(res.pop(1).(fields{j})) == 7
                    vI = volI;
                else
                    vI = volI_nn;
                end
            end
            
            r.(fields{j})(i,:) = nan(1,length(vI)); 
            % if not empty, fill data
            if sz2 > 0
                r.(fields{j})(i,vI) = res.pop(i).(fields{j});
            end
            
        end
        
        % pre compute PC for desired fields
        pcfields = {'beh_rate_adj','critp_adj','svm_rate_adj','mean_critp_adj'};
        for j = 1:length(pcfields)
            fn = sprintf('%s_PC',pcfields{j});
            r.(fn) = squeeze(PDtoPC(r.(pcfields{j})(:,[2:end]), ...
                                    r.(pcfields{j})(:,1)));
            
        end
        
        % fit logistic function to each session for each field
        logfields = {'beh_rate_adj_PC','critp_adj_PC','svm_rate_adj_PC',...
                     'mean_critp_adj_PC','auc','mean_auc'};
        for j = 1:length(logfields)
            
            fn = sprintf('%s_fit',logfields{j});
            x = r.vols_nn(i,:);
            y = r.(logfields{j})(i,:);
            xn = x(~isnan(y));
            yn = y(~isnan(y));
            lb = [min(x) .001 0 0];
            ub = [max(x) 10 1 1];
            if ~isempty(yn)
                mxslope = max(diff(yn)./diff(xn));
                [prms,mdl,thresh,sense,~,~,thresh75] = ...
                    fitLogGrid(xn,yn,[],[],[],.75,...
                               lb,ub);
                
                r.(fn).params(i,:) = prms;
                r.(fn).threshold(i,1) = thresh;
                r.(fn).sensitivity(i,1) = sense;
                r.(fn).threshold_75(i,1) = thresh75;
                r.(fn).max_slope(i,1) = mxslope;
            else
                r.(fn).params(i,:) = nan(1,4);
                r.(fn).threshold(i,1) = nan;
                r.(fn).sensitivity(i,1) = nan;
                r.(fn).threshold_75(i,1) = nan;
                r.(fn).max_slope(i,1) = nan;
            end
                
        end
            
        toc;        
        
        
    end
    
    % add in other useful info
    r.sessionID = vertcat(res.pop.sessID);
    for i = 1:length(res.pop)
        
        % mouse and contrast
        sI = ismember([sessionData.session.sessionID],r.sessionID(i));
        r.mouse{i,1} = sessionData.session(sI).mouse;
        r.contrastI(i,1) = strcmp(sessionData.session(sI).cond,'lohi');
        
    end

    save(fullfile(ops.resDir,resFile),'res','ops','r','-v7.3');
    
else
    load(fullfile(ops.resDir,resFile));

end
