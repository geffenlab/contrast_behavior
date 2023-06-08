function plot_session_summary_psych(sessID,res,spikeData,sessionData,ops)
    

clf;

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
col = s.plot.contrastColor(end,:);

%% neurogram
ax = subplot(nrows,ncols,[1 2]);
frmat = cat(1,cellRes.fr);
frmat = frmat(:,2:end) - frmat(:,1);
[~,cellSort] = sort(frmat(:,end),'descend');
frmat = frmat(cellSort,:);
cmap = zeroCMap(frmat(:),0,col,[1 1 1],[0 1 1]);
imagesc(snr(2:end),1:size(frmat,1),frmat);
h = colorbar; colormap(ax,cmap); ylabel(h,'FR_{target}-FR_{noise}');
set(gca,'xtick',snr(2:end)); set(gca,'ytick',1:size(frmat,1));
ylabel('Cell'); xlabel('Target Volume (dB SNR)');
plotPrefs; title(sprintf('%s-%d',s.mouse,s.sessionID));

%% cell waveforms
% pull out waveforms for the neurons included in these results
resCellIDs = {res.single_cell.cellID};
includedCellIDs = resCellIDs(cellI);
waveforms = spikeData.waveform.peakwaveform(...
    ismember({spikeData.cellinfo{:,7}},includedCellIDs),:);

subplot(nrows,ncols,3);
imagesc(waveforms(cellSort,:)*.195); colorbar;
title('Cell Waveforms'); xlabel('Sample Number');


%% population performance plots
subplot(nrows,ncols,4); hold on; clear ph;
xl = [snr(2) snr(end)]+[-mean(diff(snr(2:end))) mean(diff(snr(2:end)))];
plot(xl,[.5 .5],'color',[.5 .5 .5]);
ph(1) = plot(snr(2:end),PDtoPC(popRes.beh_rate_adj(2:end),popRes.beh_rate_adj(1)),...
             'color',col,'linewidth',1.5);
ph(2) = plot(snr(2:end),PDtoPC(popRes.svm_rate_adj(2:end),popRes.svm_rate_adj(1)),...
             '--','color',col,'linewidth',1);
ph(3) = plot(snr(2:end),popRes.auc,':','color',col,'linewidth',1);
ph(4) = plot(snr(2:end),PDtoPC(popRes.critp_adj(2:end),popRes.critp_adj(1)),...
             '-.','color',col,'linewidth',1);
axis tight; ylim([.4 1]);
legend(ph,'behavior','svm','pop_{auc}','pop_{crit}','location','northeastoutside');
set(gca,'xtick',snr(2:end)); set(gca,'xticklabels',num2str(snr(2:end)'));
ylabel('Percent Correct');
title('Population Results'); plotPrefs;


%% single cell performance
subplot(nrows,ncols,7:8); hold on; clear ph;
xl = [snr(2) snr(end)]+[-mean(diff(snr(2:end))) mean(diff(snr(2:end)))];
plot(xl,[.5 .5],'color',[.5 .5 .5]);
ph(1) = plot(snr(2:end),PDtoPC(popRes.beh_rate_adj(2:end),popRes.beh_rate_adj(1)),...
             'color',col,'linewidth',1.5);

if ~isempty(popRes.mean_auc)
    ph(2) = plot(snr(2:end),popRes.mean_auc,...
                 '--','color',col,'linewidth',1);
    ph(3) = plot(snr(2:end),PDtoPC(popRes.mean_critp_adj(2:end),popRes.mean_critp_adj(1)),...
                 ':','color',col,'linewidth',1);
else
    ph(2) = plot(.5,.5);
    ph(3) = plot(.5,.5);
end
axis tight; ylim([.4 1]);
legend(ph,'behavior','cell_{auc}','cell_{crit}','location','northeastoutside');
set(gca,'xtick',snr(2:end)); set(gca,'xticklabels',num2str(snr(2:end)'));
xlabel('Target Volume (dB SNR)'); ylabel('Percent Correct');
title('Target-Selective Neuron Averages'); plotPrefs;


%% session behavior metrics
subplot(nrows,ncols,6); hold on;
yyaxis left;
plot(cumsum(b.goodTrials));
ylabel('cumsum( Good Trials )');
yyaxis right;
xl = [1 length(b.goodTrials)];
plot(movmean(b.response == (b.trialType(:,1)>0),20),'k');
plot(xl,[popRes.beh_rate_adj(1),popRes.beh_rate_adj(1)],'m--');
plot(xl,[0 0],'w');
xlim(xl); ylim([0 1]);
set(gca,'ycolor','k') 

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






