function f = plot_waveforms(res);

kcols = [0 0 0; 1 0 1];
kcv = kcols(res.type,:);

subplot(2,4,[1 5]); hold on;
[rsort,sorti] = sort(res.r,'ascend');
psort = res.p(sorti);
imagesc(res.nwaveforms(sorti,:)); axis tight;
plot(xlim,[res.pcuti res.pcuti],'r','LineWidth',1.5);
plotPrefs;

subplot(2,4,[2 6]); hold on;
line(psort,1:length(psort),'Color','k','LineWidth',1.5)
plot(xlim,[res.pcuti res.pcuti],'r','LineWidth',1.5);
axis tight; plotPrefs;
ax1 = gca; 
ax1.XColor = 'r';
xlabel('p-value');
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none'); hold on;
line(rsort,1:length(rsort),'Color','b','LineWidth',1.5);
plot([0 0],[1 length(rsort)],'k:')
plot([-1 1],[res.r0 res.r0],'k:');
axis tight; plotPrefs;
ax2.XColor = 'b';
xlabel('correlation');

subplot(2,4,[4 8]);
idx = zeros(length(res.nwaveforms),1);
idx(res.include) = res.type;
[~,sorti] = sort(idx);
imagesc(res.nwaveforms(sorti,:)); set(gca,'ydir','normal');
plotPrefs;

subplot(2,4,3);
scatter(res.X(:,1)*1000,res.X(:,2)*1000,15,kcv,'.');
xlabel('Peak-Trough (ms)')
ylabel('FWHM (ms)');
plotPrefs;

subplot(2,4,7); hold on;
patchErrorBars(res.t0*1000,res.nwaveforms(idx == 0,:),[.85 .85 .85],'std');
patchErrorBars(res.t0*1000,res.nwaveforms(idx == 1,:),kcols(1,:),'std');
patchErrorBars(res.t0*1000,res.nwaveforms(idx == 2,:),kcols(2,:),'std');
xlabel('Time (ms)');
ylabel('Amplutude (au)');
plotPrefs; axis tight;