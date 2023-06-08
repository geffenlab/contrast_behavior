function plot_session_summary_offset(sessID,res,spikeData,sessionData,ops)

nrows = 3;
ncols = 3;

% get session data
sI = [sessionData.session.sessionID] == sessID;
s = sessionData.session(sI);
e = sessionData.events(sI);
b = sessionData.behavior(sI);

% get pop results for this session
popRes = res.pop([res.pop.sessID] == s.sessionID);

% get cell IDs and index into single cell results
cellIDs = {spikeData.cellinfo{[spikeData.cellinfo{:,2}] == sessID,7}};
cellI = ismember({res.single_cell.cellID},cellIDs);
cellRes = res.single_cell(cellI);

snr = popRes.snr;
offs = popRes.offs;
col = s.plot.contrastColor(end,:);

%% plot neurograms
frmat = cat(3,cellRes.fr);
frmat = frmat(2:3,:,:) - frmat(1,:,:);
[~,cellSort] = sort(squeeze(mean(frmat(2,end,:),1)),'descend');
for i = 1:2
    ax(i) = subplot(nrows,ncols,i);
    ngram = squeeze(frmat(i,:,:))';
    cmap = zeroCMap(ngram,0,s.plot.colors(i+1,:),[1 1 1],[0 1 1]);
    imagesc(1:length(offs),1:size(frmat,3),ngram(cellSort,:));
    h = colorbar; colormap(ax(i),cmap); ylabel(h,['FR_{target}-''FR_{noise}']);
    if i == 1; ylabel(sprintf('Cells\n(sorted by max FR)')); end
    set(gca,'xtick',1:length(offs)); set(gca,'xticklabels',num2str(offs'*1000));
    title(sprintf('Target Volume: %4.2f dB SNR',popRes.snr(i+1)));
    xlabel('Target Offset (ms)');
end


%% waveforms
% pull out waveforms for the neurons included in these results
resCellIDs = {res.single_cell.cellID};
includedCellIDs = resCellIDs(cellI);
waveforms = spikeData.waveform.peakwaveform(...
    ismember({spikeData.cellinfo{:,7}},includedCellIDs),:);

subplot(nrows,ncols,3);
imagesc(waveforms(cellSort,:)*.195); colorbar;
title(sprintf('%s-%d\nCell Waveforms',s.mouse,s.sessionID)); xlabel('Sample Number');


%% population performance plots
subplot(nrows,ncols,4:5); hold on; clear ph;
plot([-.025 1.025],[.5 .5],'color',[.5 .5 .5]);
plot([0 0],[.2 1],'k');
for i = 1:size(popRes.beh_rate,1)-1
    ph(1) = plot(offs,PDtoPC(popRes.beh_rate_adj(i+1,:),popRes.beh_rate_adj(1,:)),...
                'color',s.plot.colors(i+1,:),'linewidth',1.5);
    ph(2) = plot(offs,PDtoPC(popRes.svm_rate_adj(i+1,:),popRes.svm_rate_adj(1,:)),...
                 '--','color',s.plot.colors(i+1,:),'linewidth',1);
    ph(3) = plot(offs,popRes.auc(i,:),':','color',s.plot.colors(i+1,:),'linewidth',1);
    ph(4) = plot(offs,PDtoPC(popRes.critp_adj(i+1,:),popRes.critp_adj(1,:)),...
                 '-.','color',s.plot.colors(i+1,:),'linewidth',1);
end
axis tight;
legend(ph,'behavior','svm','pop_{auc}','pop_{crit}','location','northeastoutside');
set(gca,'xtick',offs); set(gca,'xticklabels',num2str(offs'*1000));
ylabel('Percent Correct');
title('Population Results'); plotPrefs;


%% single cell performance
subplot(nrows,ncols,7:8); hold on; clear ph;
plot([-.025 1.025],[.5 .5],'color',[.5 .5 .5]);
plot([0 0],[.2 1],'k');
for i = 1:size(popRes.beh_rate,1)-1
    ph(1) = plot(offs,PDtoPC(popRes.beh_rate_adj(i+1,:),popRes.beh_rate_adj(1,:)),...
                'color',s.plot.colors(i+1,:),'linewidth',1.5);
    if ~isempty(popRes.mean_auc)
        ph(2) = plot(offs,popRes.mean_auc(i,:),':','color',s.plot.colors(i+1,:),'linewidth',1);
        ph(3) = plot(offs,PDtoPC(popRes.mean_critp_adj(i+1,:),popRes.mean_critp_adj(1,:)),...
                     '-.','color',s.plot.colors(i+1,:),'linewidth',1);
    else
        ph(2) = plot(.5,.5);
        ph(3) = plot(.5,.5);
    end
end
axis tight;
legend(ph,'behavior','cell_{auc}','cell_{crit}','location','northeastoutside');
set(gca,'xtick',offs); set(gca,'xticklabels',num2str(offs'*1000));
xlabel('Target Offset (ms rel. switch)'); ylabel('Percent Correct');
title('Target-Selective Neuron Averages'); plotPrefs;


%% session behavior metrics
subplot(nrows,ncols,6); hold on;
yyaxis left;
plot(cumsum(b.goodTrials));
ylabel('cumsum( Good Trials )');
yyaxis right;
plot(movmean(b.response == (b.trialType(:,1)>0),20),'k');
ylabel('Percent Correct')
xlabel('Trial Number'); axis tight; plotPrefs;
title('Behavioral Performance');

%% plot lick RT and spike response timing
% make mean psth
cells = find(cellI);
for i = 1:length(cells)
    tpsth(i,:) = mean(res.single_cell(cells(i)).PSTH(b.trialType(:,1)>0,:));
end
mpsth = mean(tpsth,1);

% lick histogram
edges = -.5:.025:1.2;
[n,x] = hist(b.RT,edges);

subplot(nrows,ncols,9); hold on;
yyaxis left;
stairs(x,n);
ylabel('licks');
yyaxis right;
plot(ops.time,mpsth,'k');
plot([0 0],ylim,'r-');
ylabel('spks/s');
xlim([-.1 1.2]); plotPrefs;
title('Target Response Times');
xlabel('Time (s from target onset)');


            