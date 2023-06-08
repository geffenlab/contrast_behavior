function res = run_adapt(spikeData,sessionData,ops);

if ~exist('./_data/res_adapt.mat','file')

    % get neurons that have good waveforms
    included_cells = ops.include;
    spikes = spikeData.spikes(included_cells);
    cellInfo = spikeData.cellinfo(included_cells,:);

    % preallocate means
    mpsth = nan(length(spikes),length(ops.time));
    varpsth = nan(length(spikes),length(ops.time));
    for i = 1:length(spikes)
        
        % get events and session info
        s = sessionData.session(cellInfo{i,2} == [sessionData.session.sessionID]);
        e = sessionData.events(cellInfo{i,2} == [sessionData.session.sessionID]);
        b = sessionData.behavior(cellInfo{i,2} == [sessionData.session.sessionID]);
        contrastI(i) = s.stimInfo.sd(2) > s.stimInfo.sd(1);
        
        % use only noise trials
        noise_trials = b.trialType(:,1) == 0;
        
        % make psth
        psth{i} = makePSTH(spikes{i},e.trialOn(noise_trials),ops.edges);
        
        % spike mean and variance
        mpsth(i,:) = mean(psth{i},1);
        vpsth(i,:) = var(psth{i},0,1);
        
        disp(i)
        
    end

    %% standardize
    % zscored over neurons for PCA (over columns)
    time_index = ops.time_index;
    %time_index = logical(ones(size(ops.time)));
    zpsth = normalize(mpsth(:,logical(time_index)),1);
    t = ops.time(time_index);

    %% tSNE
    % activity to use for tsne
    tpsth = zpsth;

    % try tSNE
    if ~exist('Y','var')
        
        % cluster each contrast
        for j = 1:2
            
            thispsth = tpsth(contrastI == (j-1),:);
            
            Y{j} = tsne(thispsth);
            
            % create probability density and segment
            [n{j},c{j}] = hist3(Y{j},[ops.tsne_bins ops.tsne_bins]);
            ns{j} = imgaussfilt(n{j}./sum(n{j}(:)),ops.tsne_smooth);
            L{j} = double(watershed(-ns{j},8));
            
            % set up bins and labels
            [XX,YY] = meshgrid(c{j}{1},c{j}{2});
            XX = XX';
            XX = XX(:);
            YY = YY';
            YY = YY(:);
            LL = L{j}';
            LL = L{j}(:);
            XX = XX(LL>0);
            YY = YY(LL>0);
            LL = LL(LL>0);
            bincents = [XX YY];
            
            % for each neuron, find label of the nearest nonzero bin
            for i = 1:length(Y{j})
                [md,mi] = min(sqrt(sum((Y{j}(i,:) - bincents).^2,2)));
                label{j}(i) = LL(mi);
            end
        
            
        end
    end
    
    
    %% PCA
    for j = 1:2
        [coeff{j},score{j},latent{j},tsquared,explained{j}] = ...
            pca(zpsth(contrastI==(j-1),:),'NumComponents',ops.kpcs);
    end
    

    % results output
    res.cellID = {cellInfo{:,7}};
    res.contrastI = contrastI;
    res.mPSTH = mpsth;
    res.vPSTH = vpsth;
    res.tsne_embedding = Y;
    res.tsne_label = label;
    res.tsne_bins = c;
    res.tsne_n = n;
    res.tsne_pd = ns;
    res.tsne_wslabels = L;
    res.pca_coeffs = coeff;
    res.pca_score = score;
    res.pca_latent = latent;
    res.pca_explained = explained;

    save('./_data/res_adapt.mat','res','ops');

else
    load('./_data/res_adapt.mat');
end





