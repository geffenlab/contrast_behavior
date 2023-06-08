function f = plot_single_cell_offset(cellID,res,spikeData,sessionData,ops)

% get results for this cellID
rI = contains({res.single_cell.cellID},cellID);

if sum(rI) == 0
    error(sprintf(...
        'ERROR in plot_single_cell_offset.m: %s does not exist in res...\n',...
        cellID));
    f = [];
else

    f = figure; clf;
    
    % extract results for this cell
    r = res.single_cell(rI);

    % get session index and unit index for whole data set
    uI = contains(spikeData.cellinfo(:,7),cellID);
    sI = spikeData.cellinfo{uI,2} == [sessionData.session.sessionID];

    % psth and trial info
    [~,spikes,trials,PSTHs] = makePSTH(spikeData.spikes{uI},...
                                       sessionData.events(sI).trialOn + 3,...
                                       ops.edges,1);
    tt = sessionData.behavior(sI).trialType;
    [sortt,sorti] = sortrows(tt);
    [~,spkSort] = ismember(trials,sorti);



    % plot options
    pp = sessionData.session(sI).plot;
    pp.colors1 = gen_colors(sessionData.session(sI).cond,.5,.5,0,3);
    nrows = 7;

    % raster
    subplot(nrows,1,[1 2]); hold on;
    cv = pp.colors1(sortt(spkSort)+1,:);
    scatter(spikes,spkSort,2,cv,'.');
    ylim([1 length(tt)]); xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs; title(cellID,'interpreter', 'none');
    ylabel('Trials');

    % mean fr
    subplot(nrows,1,3); hold on;
    uv = unique(tt(:,1));
    for i = 1:length(uv)
        plot(ops.time,mean(PSTHs(tt(:,1)==uv(i),:),1),...
             'color',pp.colors1(i,:),'LineWidth',1);
        if i > 0
            for j = 1:length(r.offs)
                ti = ops.time > r.offs(j) & ops.time < r.offs(j)+.1;
                plot(ops.time(ti),mean(PSTHs(tt(:,1)==uv(i),ti),1),...
                     'color',pp.colors1(i,:),'LineWidth',2);
            end
        end
    end
    xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('spks/s'); xlabel('Time (s, rel. trial onset)');


    [~,spikes,trials,PSTHs] = makePSTH(spikeData.spikes{uI},...
                                       sessionData.events(sI).targOn,...
                                       ops.edges,1);

    % target triggered psth
    subplot(nrows,1,4); hold on;
    uv = unique(tt(:,1));
    for i = 1:length(uv)
        plot(ops.time,mean(PSTHs(tt(:,1)==uv(i),:),1),...
             'color',pp.colors1(i,:),'LineWidth',1);
    end
    xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('spks/s'); xlabel('Time (s, rel. target)');

    % rt histogram
    subplot(nrows,1,5); hold on;
    [n,x] = hist(sessionData.behavior(sI).RT,ops.edges);
    stairs(x,n,'color',pp.colors1(end,:),'LineWidth',1);
    xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('Lick Count'); xlabel('RT (s, rel. target)');

    % offset auc
    subplot(nrows,1,[6 7]); hold on;
    plot([r.offs(1) r.offs(end)]+[-.1 .1],[.5 .5],'k');
    for i = 1:length(r.vols)-1
        eh(i) = errorBars(r.offs+(i-1)*.01,r.auc(i,:),pp.colors1(i+1,:),...
                         [],[],squeeze(r.auc_pct(i,:,:))',[],...
                          'LineWidth',1,'Marker','.','MarkerSize',20);
    end
    plot([0 0],ylim,'k','linewidth',1);
    legend(eh,num2str(r.vols(2:end)'),'location','northeastoutside');
    set(gca,'xtick',r.offs);
    set(gca,'xticklabels',r.offs*1000);
    axis tight; plotPrefs;
    xlabel('Target Offset (ms rel. switch)');
    ylabel('AUC');

end


