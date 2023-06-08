function f = plot_single_cell_psych(cellID,res,spikeData,sessionData,ops)

% index this cell id in the res struct
rI = contains({res.single_cell.cellID},cellID);

if sum(rI) == 0
    error(sprintf(...
        'ERROR in plot_single_cell_psych.m: %s does not exist in res...\n',...
        cellID));
    f = [];
else

    clf;
    
    % extract res
    r = res.single_cell(rI);
    
    % index all units and sessions
    uI = contains(spikeData.cellinfo(:,7),cellID);
    sI = spikeData.cellinfo{uI,2} == [sessionData.session.sessionID];

    % make trial psth
    ops.edges = [-.01:.001:.1];
    ops.time = ops.edges(1:end-1) + mean(diff(ops.edges))/2;
    [PSTH,spikes,trials,PSTHs] = makePSTH(spikeData.spikes{uI},...
                                       sessionData.events(sI).targOn,...
                                       ops.edges,2);
    tt = sessionData.behavior(sI).trialType;
    [sortt,sorti] = sortrows(tt);
    [spkSortTrials,spkSort] = ismember(trials,sorti);

    % plot stuff
    pp = sessionData.session(sI).plot;
    pp.colors1 = gen_colors(sessionData.session(sI).cond,.5,.8,0);
    nrows = 6;
    ncols = 2;
    uv = unique(tt(:,1));
    xl = [-.01 .1];

    % plot target raster

    subplot(nrows,ncols,[1 2 3 4]); hold on;
    xlim(xl);
    %cv = pp.colors1(sortt(spkSort)+1,:);
    %scatter(spikes,spkSort,15,cv,'.');
    for i = 1:length(uv)
        I = [find(sortt(:,1) == (i-1),1,'first') ...
            find(sortt(:,1) == (i-1),1,'last')];
        plot(xl,[I(2) I(2)],'k');
        xx = [xl(1) xl(1) 0 0];
        yy = [I(1) I(2) I(2) I(1)];
        ph(i) = patch(xx,yy,1,'FaceColor',pp.colors1(i,:),'FaceAlpha',.5);
    end
    legend(ph,num2str(sessionData.session(sI).SNR'))
    scatter(spikes,spkSort,20,'k.');
    ylim([1 length(tt)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('Trials'); title(cellID,'interpreter','none');

    % plot mean psth
    subplot(nrows,ncols,[5 6]); hold on;
    for i = 1:length(uv)
        plot(ops.time,mean(PSTHs(tt(:,1)==uv(i),:),1),...
             'color',pp.colors1(i,:),'LineWidth',1);
    end
    xlim(xl);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('spks/s');

    % plot rt histogram
    subplot(nrows,ncols,[7 8]); hold on;
    [n,x] = hist(sessionData.behavior(sI).RT,ops.edges);
    stairs(x,n,'color',pp.contrastColor(end,:),'LineWidth',1);
    xlim(xl);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('Lick Count');
    
    % plot auc
    subplot(nrows,ncols,[9 11]); hold on;
    plot([r.vols(2) r.vols(end)]+[-2 2],[.5 .5],'k');
    plot(r.vols(2:end),r.auc,'color',pp.colors1(end,:),'linewidth',1);
    for i = 1:length(r.vols)-1
        plot(r.vols(i+1),r.auc(i),'.','markersize',20,'color',pp.colors1(i+1,:));
        errorbar(r.vols(i+1),r.auc(i),...
                 r.auc_pct(i,1)-r.auc(i),...
                 r.auc_pct(i,2)-r.auc(i),...
                 'color',pp.colors1(i+1,:),...
                 'linewidth',1)
    end
    %errorBars(r.vols(2:end),r.auc,pp.colors(end,:),[],[],r.auc_pct,[],...
%                      'LineWidth',1,'Marker','.','MarkerSize',20);
    axis tight; ylim([.4 1]); plotPrefs; axis square;
    xlabel('Target Volume (dB SNR)');
    ylabel('AUC');
    
    % plot waveform
    subplot(nrows,ncols,[10 12]); hold on;
    plot((1:63)/30e3*1000,spikeData.waveform.peakwaveform(uI,:)*.195,'k', ...
         'linewidth',1);
    xlabel('Time (ms)'); ylabel('\muV'); plotPrefs; axis square;
    
end


