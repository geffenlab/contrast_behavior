function plot_session_summary_psych_fig(sessID,res,spikeData,sessionData,ops,r,fields2use)

if ~iscell(fields2use)
    fields2use = {fields2use};
end

nrows = 3;
ncols = 3;
ms = {'o','x','*','^','d','s'};
ls = {'--',':','-.','--',':','-.'};

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
s1 = subplot(nrows,ncols,[1]);
frmat = cat(1,cellRes.fr);
frmat = frmat(:,2:end) - frmat(:,1);
[~,cellSort] = sort(frmat(:,end),'descend');
frmat = frmat(cellSort,:);
cmap = zeroCMap(frmat(:),0,col,[1 1 1],[0 1 1]);
imagesc(snr(2:end),1:size(frmat,1),frmat);
h = colorbar; colormap(s1,cmap); ylabel(h,'FR_{target}-FR_{noise}');
set(gca,'xtick',snr(2:end)); set(gca,'ytick',1:size(frmat,1));
ylabel('Cell'); xlabel('Target Volume (dB SNR)');
plotPrefs; title(sprintf('%s-%d',s.mouse,s.sessionID));

%% cell waveforms
% pull out waveforms for the neurons included in these results
resCellIDs = {res.single_cell.cellID};
includedCellIDs = resCellIDs(cellI);
waveforms = spikeData.waveform.peakwaveform(...
    ismember({spikeData.cellinfo{:,7}},includedCellIDs),:);

s2 = subplot(nrows,ncols,3);
imagesc(waveforms(cellSort,:)*.195); colorbar;
title('Cell Waveforms'); xlabel('Sample Number');


%% population performance plots
x = snr(2:end);
xl = [x(1) x(end)]+[-mean(diff(x)) mean(diff(x))]/2;
xf = linspace(snr(2),snr(end),100);

s3 = subplot(nrows,ncols,[4 5]); hold on; clear ph;
% behavior
y = PDtoPC(popRes.beh_rate_adj(2:end),popRes.beh_rate_adj(1));
[prms,mdl,thresh] = fitLogGrid(x,y);
plot(xl,[.5 .5],'k:'); %,'color',[.5 .5 .5]);
ph(1) = plot(x,y,'o','color',col,'markerfacecolor',col);
plot(xf,mdl(prms,xf),'-','color',col,'linewidth',1)
plot([thresh thresh],[.5 mdl(prms,thresh)],'-','color',col);
for i = 1:length(fields2use)
    y = PDtoPC(popRes.(fields2use{i})(2:end),popRes.(fields2use{i})(1));
    [prms,mdl,thresh] = fitLogGrid(x,y);
    ph(i+1) = plot(x+.5,y,ms{i},'color',col);
    plot(xf,mdl(prms,xf),ls{i},'color',col,'linewidth',1)
    plot([thresh thresh],[.5 mdl(prms,thresh)],ls{i},'color',col);
end
axis tight; ylim([.4 1]);
lh = [{'behavior'}, fields2use];
legend(ph,lh,'location','northeastoutside','interpreter','none');
set(gca,'xtick',x); set(gca,'xticklabels',num2str(x'));
ylabel('Percent Correct'); xlabel('Target Volume (dB SNR)');
plotPrefs;

% histogram for top performance
normp = popRes.projection;
include = ~b.abort & b.goodTrials;
vols = b.trialType(include,1);
[pred, crit, rule] = crit_classifier(normp,vols);

s4 = subplot(nrows,ncols,6); hold on;
edges = linspace(0,1,20);
histogram(normp(vols==0),edges,'facecolor',[.5 .5 .5],'normalization','probability');
histogram(normp(vols==6),edges,'facecolor',col,'normalization','probability');
plot([crit crit],ylim,'k','linewidth',1);
xlabel('CD Projection Value'); axis tight;
ylabel('p(Proj. Val.)'); xlim([edges(1) edges(end)]);
plotPrefs;



%% session behavior metrics
c2 = [28,144,153] / 256;
s5 = subplot(nrows,ncols,7:8); hold on;
xl = [1 length(b.goodTrials)];
plot(xl,[.5 .5],'k:'); %,'color',[.5 .5 .5]);
plot(cumsum(b.goodTrials)./sum(b.goodTrials),'color',c2,'LineWidth',1);
plot(movmean(b.response == (b.trialType(:,1)>0),20),'k','LineWidth',1);
plot(xl,[popRes.beh_rate_adj(1),popRes.beh_rate_adj(1)],'m--','LineWidth',1);
ylabel('Probability/CDF'); ylim([0 1]);
xlabel('Trial Number'); xlim(xl); plotPrefs;
legend('chance','cdf(Good Trials)','Percent Correct','Mean False Alarm Rate',...
       'location','northeastoutside');


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

s6 = subplot(nrows,ncols,9); hold on;
yyaxis left;
plot(ops.time,mpsth,'k','LineWidth',1);
ylabel('spks/s');
yyaxis right;
plot([0 0],[0 max(n)].*1.1,'r-','LineWidth',1);
stairs(x,n,'color',c2,'LineWidth',1);
ylabel('licks'); axis tight;
xlim([-.1 1.2]); plotPrefs;
title('Target Response Times');
xlabel('Time (s from target onset)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = c2;

s3.Position(3) = s3.Position(3) * .66;
s1.Position(3) = s3.Position(3);
s5.Position(3) = s6.Position(3);








