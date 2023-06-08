function [res,data,ops,exitflag] = lnmodelcv(data,res,c,plotDir,ops,debugMode)

% this function runs an LN model on spike data, fitting an STA and
% nonlinearity to the provided cell during low and high contrast
% periods

% c = 22
% c = 56;
% c = 3
% c = 6
% c = 44
% c = 34

% in debug mode, load simulated neuron
if exist('debugMode','var') & debugMode
    fprintf('\n***\nWARNING: RUNNING IN DEBUG MODE!!\n***\n');
    dataDir = './_data';
    fn = 'simulated_neuron.mat';
    load(fullfile(dataDir,fn));
    c = 1;
end

exitflag = false;

% models to run
models = ops.models;

% current stimulus file
fp = strsplit(data.behavior.stimInfo.stim,'\');
fn = fp{end};

% check if the stimulus is appropriate for this LN analysis
fileList = dir(fullfile(pwd,'_stimuli'));
if any(contains({fileList.name},fn))
    
    % save original ops in a separate structure for reference
    ln = ops;

    if length(res) >= c & isfield(res(c),'noiseRatio')
        % if this analysis is already run
        fprintf('\tAlready ran analysis: skipping...\n');
        exitflag = true;
        return;
                
    else
        
        fprintf('\t\tLN MODEL\n');
                
        
        %% DATA FORMATTING
        % load and process the stimuli if not done already            
        formatStimAndSpikes;
        
        % matched stim index
        I = data.stim.index(7,:) == 1;
        
        % check if there aren't many spikes
        if sum(spikes(I)) < 10
            fprintf('\tLess than 10 spikes, skipping... \n');
            exitflag = true;
            return;
        end
  
  
        
        %% CROSSVAL SETUP
        % counterbalancing
        cvI = crossvalind('Kfold',data.behavior.trialType(:,3), ...
                          ops.kfolds);
        
        if any(contains(models,'glm'))
            % crossvalidate GLM regularization parameter (lambda) for
            % all matched data
            fprintf('\t\tFinding glmnet lambda... '); tic;
            X = makeDesignMatrix(data.stim.db(:,I),length(ops.t));
            CVerr = cvglmnet(X',spikes(I)','gaussian',[],[],10);
            lam = CVerr.lambda_min; % lambda value to use
            toc;
        end
        
        % preallocate and clear
        for i = 1:length(models)
            res(c).(models{i}).y = nan(length(data.stim.index),1);
            res(c).(models{i}).yhat = nan(length(data.stim.index),1);
            res(c).(models{i}).ylin = nan(length(data.stim.index),1);
            static_lbl = sprintf('%s_static',models{i});
            res(c).(static_lbl).y = nan(length(data.stim.index),1);
            res(c).(static_lbl).yhat = nan(length(data.stim.index),1);
            res(c).(static_lbl).ylin = nan(length(data.stim.index),1);
                        
        end
        
        
        
        %% CROSSVAL LOOP
        fprintf('\t\tCounterbalanced fit: cv fold '); tic;
        for i = 1:ops.kfolds
            
            fprintf('%d ',i);
            
            % train and test index
            trainI = ismember(data.stim.index(1,:),find(cvI ~= i));
            testI = ismember(data.stim.index(1,:),find(cvI == i));
            
            
            %% FILTERS
            % make design matrix from matched data in each contrast
            I = data.stim.index(7,:) == 1 & trainI;  
            X = makeDesignMatrix(data.stim.db(:,I),length(ops.t));
            
            if any(contains(models,'sta'))
                % spike triggered average
                sta = genSTA(find(spikes>0 & I),data.stim.db,ln.w,ln.fs);
                sta = sta - mean(sta,2); % subtract mean over time
                filter.sta = normSTA(sta);
            end
            
            if any(contains(models,'rc'))
                % RC-STA (not normalized)
                rc = rcSTA(spikes(I),X,length(ops.t),false);
                rc = rc - mean(rc,2); % subtract mean over time
                filter.rc = normSTA(rc);
            end
            
            if any(contains(models,'nrc'))
                % NRC-STA (normalized)
                nrc = rcSTA(spikes(I),X,length(ops.t),true);
                nrc = nrc - mean(nrc,2); % subtract mean over time
                filter.nrc = normSTA(nrc);
            end
            
            if any(contains(models,'glm'))
                % GLM STRF w. regularization
                options.lambda = lam;
                %options.alpha = 0;
                fit = glmnet(X',spikes(I)','gaussian',options);
                glm = fliplr(reshape(fit.beta,length(ops.f),[]));
                filter.glm = normSTA(glm);
            end
            
            
            
            %% static nonlinearity
            trainops = ops;
            trainops.index = [trainI & ...
                              data.stim.index(7,:) == 1];
            trainops.sigma = [];
            
            % testing
            testops = trainops;
            testops.index = [testI & ...
                             data.stim.index(7,:) == 1];
            testops.sigma = ops.sigma;
            
            % train and test each model
            for k = 1:length(models)
                
                m = models{k};
                
                %% TRAIN
                trainops.sta = filter.(m);
                train.(m) = fitLN(data.stim.db,spikes,trainops);

                %% TEST
                test.(m) = predictLN(data.stim.db,spikes,train.(m),testops);
                
                %% RES
                lbl = sprintf('%s_static',m);
                res(c).(lbl).xy{i} = test.(m).xy;
                res(c).(lbl).my(i,:) = [min(test.(m).xy(:,1)) max(test.(m).xy(:,2))];
                res(c).(lbl).ahat(i,:) = test.(m).ahat;
                res(c).(lbl).r(i) = test.(m).r;
                res(c).(lbl).mse(i) = test.(m).MSE;
                res(c).(lbl).filter(:,:,i) = test.(m).sta;
                res(c).(lbl).y(testops.index) = test.(m).y;
                res(c).(lbl).yhat(testops.index) = test.(m).yhat;
                res(c).(lbl).ylin(testops.index) = test.(m).ylintest;
                res(c).(lbl).mdl = test.(m).model;
                
            end             
            
            
            %% CONTRAST LOOP FOR GC MODEL
            % for each contrast
            for j = 1:2
                
                %% training and testing ops
                % training
                trainops = ops;
                trainops.index = [trainI & ...
                                  data.stim.index(2,:) == j & ...
                                  data.stim.index(7,:) == 1];
                trainops.sigma = [];
                
                % testing
                testops = trainops;
                testops.index = [testI & ...
                                 data.stim.index(2,:) == j & ...
                                 data.stim.index(7,:) == 1];
                testops.sigma = ops.sigma;
                
                % train and test each model
                for k = 1:length(models)
                    
                    m = models{k};
                    
                    %% TRAIN
                    trainops.sta = filter.(m);
                    train.(m) = fitLN(data.stim.db,spikes,trainops);

                    %% TEST
                    test.(m) = predictLN(data.stim.db,spikes,train.(m),testops);
                    
                    %% RES
                    res(c).(m).xy{i,j} = test.(m).xy;
                    res(c).(m).my(i,j,:) = [min(test.(m).xy(:,1)) max(test.(m).xy(:,2))];
                    res(c).(m).ahat(i,j,:) = test.(m).ahat;
                    res(c).(m).r(i,j) = test.(m).r;
                    res(c).(m).mse(i,j) = test.(m).MSE;
                    res(c).(m).filter(:,:,i) = test.(m).sta;
                    res(c).(m).y(testops.index) = test.(m).y;
                    res(c).(m).yhat(testops.index) = test.(m).yhat;
                    res(c).(m).ylin(testops.index) = test.(m).ylintest;
                    res(c).(m).mdl = test.(m).model;
                    
                end
                
            end
        end
        toc;
                
        
        %% RESHAPE BY TRIAL    
        % find max trial length to use for padding
        uTrials = unique(data.stim.index(1,:));
        for i = 1:length(uTrials)
            n(i) = sum(data.stim.index(1,:)==i);
        end
        nMax = max(n);
        nMin = min(n);
                
        % preallocate
        allocate = nan(length(uTrials),nMax);
        
        % get values for each model
        labels = models;
        labels{end+1} = sprintf('%s_static',models{1});
        labels{end+1} = sprintf('%s_static',models{2});
        labels{end+1} = sprintf('%s_static',models{3});
        for i = 1:length(labels)  
            
            m = labels{i};

            res(c).(m).yshape = allocate;
            res(c).(m).yhatshape = allocate;
            res(c).(m).ylinshape = allocate;
            
            for j = 1:length(uTrials)
                
                I = find(data.stim.index(1,:)==j);
                n = length(I);
                
                % add to matrix
                res(c).(m).yshape(j,1:n) = res(c).(m).y(I);
                res(c).(m).ylinshape(j,1:n) = res(c).(m).ylin(I);
                res(c).(m).yhatshape(j,1:n) = res(c).(m).yhat(I);

            end
                        
        end
        
        % time index
        timeI = ops.bin/2:ops.bin:nMax*ops.bin-(ops.bin/2);
        contrastI = timeI > 3;
        
        
        
        
        %% CORRELATIONS AND NOISE RATIO
        % mean correlations
        for i = 1:length(labels)
            for j = 1:2
                I = contrastI+1 == j;
                y = nanmean(res(c).(labels{i}).yshape(:,I))';
                yhat = nanmean(res(c).(labels{i}).yhatshape(:,I))';
                
                res(c).(labels{i}).meanr(j) = corr(y,yhat,'rows','complete');
            end
            res(c).(labels{i}).allr = corr(nanmean(res(c).(labels{i}).yshape)',...
                                         nanmean(res(c).(labels{i}).yhatshape)',...
                                           'rows','complete');
            res(c).(labels{i}).allmse = nanmean((nanmean(res(c).(labels{i}).yshape) - ...
                                           nanmean(res(c).(labels{i}).yhatshape)).^2);
        end
            

        % make a raster
        trialEdges = -.1:.005:5.1;
        trialTime = trialEdges(1:end-1) + .005/2;
        [noiseSortT,noiseSortI] = ...
            sortrows([data.behavior.trialType(:,3) data.behavior.trialType(:,1)]);
        [trialPSTH,trialRaster,trialTrials,trialPSTHsm] = ...
            makePSTH(data.spikes{c},data.events.trialOn,trialEdges,5);
        [~,noiseSortTrials] = ismember(trialTrials,noiseSortI);
        
        
        % compute noise power, excluding transition periods and
        % excluding target trials
        I = data.behavior.trialType(:,1) == 0;
        for i = 1:2
            NR(:,i) = responsePower(trialPSTH(I,contrastI+1 == i), ...
                                  data.behavior.trialType(I,3));
        end
        res(c).noiseRatio = NR;
        
        % add some other stuff
        res(c).model = test.glm.model;
        res(c).ops = ops;
        res(c).models = models;
        res(c).modelStrs = ops.modelStrings;
        res(c).maxTrialSamps = nMax;

        
        
        
        
        %% PLOTTING
        if exist('plotDir','var') & ~isempty(plotDir)
                        
            tStr = sprintf('%s Cell %d (%s) | %s-%s',data.params.plot.titleString,...
                           c,data.cellInfo{c,6},data.params.cond, ...
                           data.params.task);
            
            if exist('debugMode','var') & debugMode
                tStr = sprintf('simulated neuron:\n p1=[%03.2f %03.2f %03.2f %03.2f], p2=[%03.2f %03.2f %03.2f %03.2f]',...
                               data.debug.p0(1,:)', data.debug.p0(2,:)');
            end
            
            % compute plot size based on models
            modelN = length(models);
            plotN = 2 + length(models);
            
            % plot raster
            subplot(4,plotN,[1 2 plotN+1 plotN+2])
            hold on
            scatter(trialRaster,noiseSortTrials,10,'.k')
            xlim([trialEdges(1) trialEdges(end)]);
            p = patch([0 0 3 3],length(noiseSortT)+[0 20 20 0],1);
            p.EdgeAlpha = 0;
            p.FaceColor = data.params.plot.contrastColor(1,:);
            p.FaceAlpha = .75;
            p = patch([3 3 data.params.offsets(end)+4 data.params.offsets(end)+4],...
                      length(noiseSortT)+[0 20 20 0],1);
            p.EdgeAlpha = 0;
            p.FaceColor = data.params.plot.contrastColor(2,:);
            p.FaceAlpha = .75;
            axis tight
            plot([0 0],ylim,'k','LineWidth',1);
            plot([3 3],ylim,'k','LineWidth',1);
            clear yticks;
            for i = 1:length(unique(noiseSortT(:,1)))
                ind = [find(noiseSortT(:,1) == i,1,'first') ...
                       find(noiseSortT(:,1) == i,1,'last')];
                yticks(i) = mean(ind);
                plot(xlim,[ind(2) ind(2)],'k');
            end
            set(gca,'ytick',yticks);
            set(gca,'yticklabels',{'FN1','FN2','FN3','FN4','FN5'});
            plotPrefs; title(sprintf('%s\nTrial-triggered Responses',tStr));
            
            %% plot model components
            % color limit relative to first model
            clim = [min(min(mean(res(c).(models{1}).filter,3))) ...
                    max(max(mean(res(c).(models{1}).filter,3)))];
            for i = 1:length(models)
                % filter plot
                s1(i) = subplot(4,plotN,2+[0 plotN]+i);
                plotSTA(ops.t,ops.f,fliplr(mean(res(c).(models{i}).filter,3)),[],clim);
                colorbar; xlabel('Time (s)'); ylabel('Frequency (kHz)');
                title(res(c).modelStrs{i});
                plotPrefs; axis square;
                
                tmp = vertcat(res(c).(models{i}).xy{:});
                xr(i,:) = [min(tmp(:,1)) max(tmp(:,1))];
                yr(i,:) = [min(tmp(:,2)) max(tmp(:,2))];
                
                % nonlinearities
                subplot(4,plotN,2+[plotN*2 plotN*3]+i); hold on;
                
                % static model
                lbl = sprintf('%s_static',models{i});
                tmp = vertcat(res(c).(lbl).xy{:});
                scatter(tmp(:,1),tmp(:,2),5,...
                        'MarkerFaceColor',[.7 .7 .7],...
                        'MarkerEdgeColor','none',...
                        'MarkerFaceAlpha',.5)
                x = sort(tmp(:,1));
                for k = 1:size(res(c).(lbl).ahat,1)
                    plot(x,res(c).model(res(c).(lbl).ahat(k,:),x),...
                         'Color',[.7 .7 .7],...
                         'LineWidth',1);  
                end
                
                % plot nonlinearity scatter and fits
                for j = 1:2
                    tmp = vertcat(res(c).(models{i}).xy{:,j});
                    s = scatter(tmp(:,1),tmp(:,2),5);
                    s.MarkerFaceColor = data.params.plot.contrastColor(j,:);
                    s.MarkerEdgeColor = 'none';
                    s.MarkerFaceAlpha = .5;
                    x = sort(tmp(:,1));
                    for k = 1:size(res(c).(models{i}).ahat,1)
                        plot(x,res(c).model(res(c).(models{i}).ahat(k,j,:),x),...
                             'Color',data.params.plot.contrastColor(j,:),...
                             'LineWidth',1);  
                    end
                    if exist('debugMode','var') & debugMode
                        plot(x,data.debug.model(data.debug.p0(j,:),x),'--', ...
                             'Color',data.params.plot.contrastColor(j,:),...
                             'LineWidth',1);
                    end
                end
                xlabel('Filter Response'); ylabel('FR (Hz)');
                plotPrefs; axis square;

                
            end
            
            xl = round(max(abs(xr(:))))+1;
            yl = round(max(abs(yr(:))))+1;
            
            for i = 1:length(models)
                s2 = subplot(4,plotN,2+[plotN*2 plotN*3]+i);
                xlim([-xl xl]); ylim([0 yl]);
                s2.Position = [s2.Position(1:2) s1(i).Position(3:4)];
            end
            
            
            
            %% plot predictions
            s = subplot(4,plotN,[1 2] + plotN*3); hold on;
            clear p;
            labels = {'glm_static','glm'};
            for i = 1:length(labels)
                p(i) = plot(timeI,nanmean(res(c).(labels{i}).yhatshape),'LineWidth',1);
            end
            p(end+1) = plot(timeI,nanmean(res(c).(models{1}).yshape),'k');
            plot([0 0],ylim,'k','LineWidth',1);
            plot([3 3],ylim,'k','LineWidth',1);
            lstr = labels; lstr{end+1} = 'observed';
            legend(p,lstr,'location','nw','interpreter','none');
            title(sprintf('Model Predictions\nr_{GC}=%03.2f, r_{static}=%03.2f',...
                          res(c).glm.allr,res(c).glm_static.allr)); 
            xlabel('Time (s)'); ylabel('FR (Hz)');
            axis tight; xlim([trialEdges(1) trialEdges(end)]);
            plotPrefs;
            
                                        
        end
        
    end
    
    
else
    fprintf(['\tAppropriate stim not found, skipping LN analysis...\' ...
             'n']);
    ops = [];
    exitflag = true;
    
end

