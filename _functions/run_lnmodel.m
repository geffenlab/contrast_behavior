function [res] = run_lnmodel(spikeData,sessionData,ops)

%% function [res] = run_lnmodel(spikeData,sessionData,ops)
resFile = ops.resFile;
resFN = fullfile(ops.resDir,resFile);

fprintf('Checking stim files... '); tic;
% make a big cell array with the stimulus for each cell
incl = find(ops.include);
for i = 1:length(incl)
    % find the session for this cell
    sI = find(vertcat(sessionData.session.sessionID) == ...
              spikeData.cellinfo{incl(i),2});
    sessID = sessionData.session(sI).stimInfo.IDsess;
    stimfile = fullfile('./_data','_spectrograms',[sessID '-spec.mat']);
    stimExist(i) = exist(stimfile);
end
stimcells = false(size(ops.include));
stimcells(incl) = stimExist>0;
toc;

% of included neurons, find ones with the right stimulus and with
% higher fr
included_cells = stimcells;
spikes = spikeData.spikes(included_cells);
cellInfo = spikeData.cellinfo(included_cells,:);
if sum(stimcells) == 0
    error('No matching stimuli found... check if ./_data/_spectrograms is empty')
end

% session for the first neuron
last_sI = 0;

% models to run
models = ops.models;

% figure handling
sz = [1500 400];
f1 = figure(1); set(f1,'visible',ops.fig_visible,'Position',[0 0 sz]);


%% single neurons
clear sI;
for c = 1:length(spikes)

    % c = 122;

    % cell info
    fprintf('CELL %d/%d... \n',c,length(spikes));
    fn = sprintf('./_data/_lnmodel/%s.mat',cellInfo{c,7});
    sI(c) = find(vertcat(sessionData.session.sessionID) == ...
                 cellInfo{c,2});

    t0 = tic;

    if ~exist(fn,'file') & ~exist(resFN, 'file')


        % format spikes and stimuli (only if cell is in a new session)
        stimInfo = sessionData.session(sI(c)).stimInfo;
        behavior = sessionData.behavior(sI(c));
        ops.f = stimInfo.freqs;
        stimInfo.offsets = sessionData.session(sI(c)).offsets;
        events = sessionData.events(sI(c)).trialOn;
        [stim,spec,y0,ops] = formatStimAndSpikes(...
            spikes{c},events,stimInfo,behavior,ops,last_sI==sI(c));


        % index for matched stimuli
        I = stim.index(7,:) == 1;
        if sum(y0(I)) < 10
            fprintf('\tLess than 10 spikes, skipping... \n');
            continue;
        end

        % design matrix
        S = stim.db - mean(stim.db,2);
        S = S - mean(S(:));
        X0 = lagDesignMatrix(S,length(ops.t));

        if (2+2) == 5
            cvfit = cvglmnet(X0',y0','gaussian',[],[],10);
            coeffs = cvglmnetCoef(cvfit,'lambda_1se');
            imagesc(reshape(coeffs(2:end),length(ops.f),[]))
        end


        %% CROSSVAL SETUP
        % counterbalancing
        cvI = crossvalind('Kfold',behavior.trialType(:,3), ...
                          ops.kfolds);

        if any(contains(models,'glm'))
            % crossvalidate GLM regularization parameter (lambda) for
            % all matched data
            fprintf('\tFinding glmnet lambda... '); tic;
            cvfit = cvglmnet(X0(:,I)',y0(I)','gaussian',[],[],10);
            coeffs = cvglmnetCoef(cvfit,'lambda_1se');
            lam = cvfit.lambda_1se; % lambda value to use
            toc;
        end

        % preallocate and clear
        clear res;
        for i = 1:length(models)
            res.(models{i}).y = nan(length(stim.index),1);
            res.(models{i}).yhat = nan(length(stim.index),1);
            res.(models{i}).ylin = nan(length(stim.index),1);
            static_lbl = sprintf('%s_static',models{i});
            res.(static_lbl).y = nan(length(stim.index),1);
            res.(static_lbl).yhat = nan(length(stim.index),1);
            res.(static_lbl).ylin = nan(length(stim.index),1);
        end

        %% CROSSVAL LOOP
        fprintf('\tCounterbalanced fit: cv fold '); tic;
        for i = 1:ops.kfolds

            fprintf('%d ',i);

            % train and test index
            trainI = ismember(stim.index(1,:),find(cvI ~= i));
            testI = ismember(stim.index(1,:),find(cvI == i));

            %% FILTERS
            I = stim.index(7,:) == 1 & trainI;
            X = makeDesignMatrix(stim.db(:,I),length(ops.t));

            if any(contains(models,'sta'))
                % spike triggered average
                sta = genSTA(find(y0>0 & I),stim.db,ops.w,ops.fs);
                sta = sta - mean(sta,2); % subtract mean over time
                filter.sta = normSTA(sta);
            end

            if any(contains(models,'rc'))
                % RC-STA (not normalized)
                rc = rcSTA(y0(I),X,length(ops.t),false);
                rc = rc - mean(rc,2); % subtract mean over time
                filter.rc = normSTA(rc);
            end

            if any(contains(models,'NRC'))
                % NRC-STA (normalized)
                nrc = rcSTA(y0(I),X,length(ops.t),true);
                nrc = nrc - mean(nrc,2); % subtract mean over time
                filter.NRC = normSTA(nrc);
            end

            if any(contains(models,'glm'))
                % GLM STRF w. regularization
                %options.alpha = 1;
                %fit = cvglmnet(X0(:,I)',y0(I)','gaussian',options,[],5);
                %coeffs = cvglmnetCoef(fit,'lambda_1se');
                %glm = fliplr(reshape(coeffs(2:end),length(ops.f),[]));

                options.lambda = lam;
                fit = glmnet(X',y0(I)','gaussian',options);
                glm = fliplr(reshape(fit.beta,length(ops.f),[]));
                filter.glm = normSTA(glm);
            end



            %% static nonlinearity
            trainops = ops;
            trainops.index = [trainI & ...
                              stim.index(7,:) == 1];
            trainops.sigma = [];

            % testing
            testops = trainops;
            testops.index = [testI & ...
                             stim.index(7,:) == 1];
            testops.sigma = ops.sigma;

            % train and test each model
            for k = 1:length(models)

                m = models{k};

                %% TRAIN
                trainops.sta = filter.(m);
                train.(m) = fitLN(stim.db,y0,trainops);

                %% TEST
                test.(m) = predictLN(stim.db,y0,train.(m),testops);

                %% RES
                lbl = sprintf('%s_static',m);
                res.(lbl).xy{i} = test.(m).xy;
                res.(lbl).my(i,:) = [min(test.(m).xy(:,1)) max(test.(m).xy(:,2))];
                res.(lbl).ahat(i,:) = test.(m).ahat;
                res.(lbl).r(i) = test.(m).r;
                res.(lbl).mse(i) = test.(m).MSE;
                res.(lbl).filter(:,:,i) = test.(m).sta;
                res.(lbl).y(testops.index) = test.(m).y;
                res.(lbl).yhat(testops.index) = test.(m).yhat;
                res.(lbl).ylin(testops.index) = test.(m).ylintest;
                res.(lbl).mdl = test.(m).model;

            end



            %% CONTRAST LOOP FOR GC MODEL
            % for each contrast
            for j = 1:2

                %% training and testing ops
                % training
                trainops = ops;
                trainops.index = [trainI & ...
                                  stim.index(2,:) == j & ...
                                  stim.index(7,:) == 1];
                trainops.sigma = [];

                % testing
                testops = trainops;
                testops.index = [testI & ...
                                 stim.index(2,:) == j & ...
                                 stim.index(7,:) == 1];
                testops.sigma = ops.sigma;

                % train and test each model
                for k = 1:length(models)

                    m = models{k};

                    %% TRAIN
                    trainops.sta = filter.(m);
                    train.(m) = fitLN(stim.db,y0,trainops);

                    %% TEST
                    test.(m) = predictLN(stim.db,y0,train.(m),testops);

                    %% RES
                    res.(m).xy{i,j} = test.(m).xy;
                    res.(m).my(i,j,:) = [min(test.(m).xy(:,1)) max(test.(m).xy(:,2))];
                    res.(m).ahat(i,j,:) = test.(m).ahat;
                    res.(m).r(i,j) = test.(m).r;
                    res.(m).mse(i,j) = test.(m).MSE;
                    res.(m).filter(:,:,i) = test.(m).sta;
                    res.(m).y(testops.index) = test.(m).y;
                    res.(m).yhat(testops.index) = test.(m).yhat;
                    res.(m).ylin(testops.index) = test.(m).ylintest;
                    res.(m).mdl = test.(m).model;

                end

            end
        end
        toc;


        %% RESHAPE BY TRIAL
        % find max trial length to use for average PSTH
        nMax = max(spec.trial_actual_length + 2*spec.pad_samps);

        % preallocate
        allocate = nan(length(behavior.trialType),nMax);

        % get values for each model
        labels = models;
        for i = 1:length(labels)
            labels{end+1} = sprintf('%s_static',models{i});
        end
        for i = 1:length(labels)

            m = labels{i};

            res.(m).yshape = allocate;
            res.(m).yhatshape = allocate;
            res.(m).ylinshape = allocate;

            for j = 1:length(behavior.trialType)

                I = find(stim.index(1,:)==j);
                n = length(I);

                % add to matrix
                res.(m).yshape(j,1:n) = res.(m).y(I);
                res.(m).ylinshape(j,1:n) = res.(m).ylin(I);
                res.(m).yhatshape(j,1:n) = res.(m).yhat(I);

            end

        end

        timeI = -ops.pad : 1/ops.fs : (nMax/ops.fs-ops.pad);
        timeI = timeI(1:end-1) + 1/ops.fs/2;
        contrastI = timeI > 3;


        %% CORRELATIONS AND NOISE RATIO
        for i = 1:length(labels)
            for j = 1:2
                I = contrastI+1 == j;
                y = nanmean(res.(labels{i}).yshape(:,I))';
                yhat = nanmean(res.(labels{i}).yhatshape(:,I))';

                res.(labels{i}).meanr(j) = corr(y,yhat,'rows','complete');
            end
            res.(labels{i}).allr = corr(nanmean(res.(labels{i}).yshape)',...
                                           nanmean(res.(labels{i}).yhatshape)',...
                                           'rows','complete');
            res.(labels{i}).allmse = nanmean((nanmean(res.(labels{i}).yshape) - ...
                                                 nanmean(res.(labels{i}).yhatshape)).^2);
        end

        % make a raster
        trialEdges = -.1:.005:5.1;
        trialTime = trialEdges(1:end-1) + .005/2;
        [noiseSortT,noiseSortI] = ...
            sortrows([behavior.trialType(:,3) behavior.trialType(:,1)]);
        [trialPSTH,trialRaster,trialTrials,trialPSTHsm] = ...
            makePSTH(spikes{c},events,trialEdges,5);
        [~,noiseSortTrials] = ismember(trialTrials,noiseSortI);


        % compute noise power, excluding transition periods and
        % excluding target trials
        I = behavior.trialType(:,1) == 0;
        for i = 1:2
            NR(:,i) = responsePower(trialPSTH(I,contrastI+1 == i), ...
                                    behavior.trialType(I,3));
        end
        res.noiseRatio = NR;

        % add some other stuff
        res.model = test.glm.model;
        res.ops = ops;
        res.models = models;
        res.modelStrs = ops.modelStrings;
        res.maxTrialSamps = nMax;
        res.cellInfo = cellInfo(c,:);
        res.cellID = cellInfo{c,7};
        res.sessionID = sessionData.session(sI(c)).sessionID;
        res.cond = sessionData.session(sI(c)).cond;
        res.sessionIndex = sI(c);
        res.spec_file = spec.spec_file;
        res.timeIndex = timeI;
        res.contrastIndex = contrastI;
        res.behavior = behavior;
        res.events = sessionData.events(sI(c));
        res.ops = ops;


        %% save
        r = res;
        save(fn,'r');


        %% plot results for this neuron
        set(0,'CurrentFigure',f1); clf(f1);
        plot_lnmodel_cell(res,spikes(c),models);
        plot_place = './_plots/_lnmodel';
        if ~exist(plot_place,'dir')
            mkdir(plot_place);
        end
        saveFigPDF(f1,sz,sprintf('%s/%s_lnmodels.pdf', plot_place, res.cellID));

    else

        fprintf('\tFile found, skipping...\n');
        %load(fn,'r');
        %res = r;

    end

    fprintf('\tRuntime: '); toc(t0);


end

if ~exist(resFN)

    %% build a results file
    models2use = {'sta','sta_static','NRC','NRC_static','glm', ...
                  'glm_static'};
    fields2use = {'meanr','allr','allmse','ahat','filter'};

    clear res
    for c = 1:numel(sI)
        
        % load
        fn = sprintf('./_data/_lnmodel/%s.mat',cellInfo{c,7});
        fprintf('Res for %s (%d/%d)... ',cellInfo{c,7},c,numel(sI)); tic;
        load(fn)
        
        % for each model
        for i = 1:length(models2use)
            % for each desired field
            for j = 1:length(fields2use)
                fieldname = sprintf('%s_%s',models2use{i},fields2use{j});
                if contains(fields2use{j},'filter')
                    res(c).(fieldname) = squeeze(mean(r.(models2use{i}).(fields2use{j}),3));
                else
                    res(c).(fieldname) = squeeze(mean(r.(models2use{i}).(fields2use{j}),1));                    
                end
            end
        end
        
        % individual stuff
        res(c).cellInfo = r.cellInfo;
        res(c).cellID = r.cellID;
        res(c).sessionID = r.sessionID;
        res(c).cond = r.cond;
        res(c).sessionIndex = r.sessionIndex;
        res(c).noiseRatio = r.noiseRatio;
        res(c).ops = r.ops;
        res(c).task = r.behavior.task;
        res(c).behavior.goodTrials = r.behavior.goodTrials;
        res(c).behavior.response = r.behavior.response;
        res(c).behavior.trialType = r.behavior.trialType;
        res(c).behavior.abort = r.behavior.abort;
        toc;
        
    end

    ops.models2use = models2use;
    ops.fields2use = fields2use;

    save(resFN,'res','ops');

else
    
    load(resFN);
    
end
