function plot_lnmodel_cell(res,spks,mdls2use)


%% setup
mdl_colors = [27,158,119;...
              217,95,2;...
              117,112,179;...
              231,41,138;...
              102,166,30;...
              230,171,2;...
              166,118,29;...
              102,102,102] / 256;
mdl_colors_lite = [102,194,165;...
                   252,141,98;...
                   141,160,203;...
                   231,138,195;...
                   166,216,84;...
                   255,217,47;...
                   229,196,148;...
                   179,179,179] / 256;
if contains(res.cond,'hilo')
    cc = [1 0 0; 0 0 1];
else
    cc = [0 0 1; 1 0 0];
end
nrows = 4;
ncols = 2 + length(mdls2use);

if iscell(spks)
    spks = spks{1};
end

%% raster
% compute
trialEdges = -.1:.005:5.1;
trialTime = trialEdges(1:end-1) + .005/2;
[noiseSortT,noiseSortI] = ...
    sortrows([res.behavior.trialType(:,3) res.behavior.trialType(:,1)]);
[trialPSTH,trialRaster,trialTrials,trialPSTHsm] = ...
    makePSTH(spks,res.events.trialOn,trialEdges,5);
[~,noiseSortTrials] = ismember(trialTrials,noiseSortI);

% plot raster
subplot(4,ncols,[1 2 ncols+1 ncols+2])
hold on
scatter(trialRaster,noiseSortTrials,10,'.k')
xlim([trialEdges(1) trialEdges(end)]);
p = patch([0 0 3 3],length(noiseSortT)+[0 20 20 0],1);
p.EdgeAlpha = 0;
p.FaceColor = cc(1,:);
p.FaceAlpha = .75;
p = patch([3 3 4.5 4.5],...
          length(noiseSortT)+[0 20 20 0],1);
p.EdgeAlpha = 0;
p.FaceColor = cc(2,:);
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
set(gca,'yticklabels',{'Scene 1','Scene 2','Scene 3','Scene 4','Scene 5'});
plotPrefs; 
title(sprintf('%s\nTrial-triggered Responses',res.cellID),...
      'interpreter','none');


%% predictions
s = subplot(4,ncols,[1 2] + ncols*3); hold on;
clear p lstr;
fn = fieldnames(res);
labels = fn(contains(fn,mdls2use));
for i = 1:length(labels)
    p(i) = plot(res.timeIndex,nanmean(res.(labels{i}).yhatshape),...
                'Color',mdl_colors(i,:));
    
    lstr{i} = sprintf('%s (r=%0.3f)',labels{i},res.(labels{i}).allr);
end
p(end+1) = plot(res.timeIndex,nanmean(res.(labels{1}).yshape),'k','linewidth',1);
plot([0 0],ylim,'k','LineWidth',1);
plot([3 3],ylim,'k','LineWidth',1);
lstr{end+1} = 'spikes';
legend(p,lstr,'location','nw','interpreter','none');
title('Model Predictions'); 
xlabel('Time (s)'); ylabel('FR (Hz)');
axis tight; xlim([trialEdges(1) trialEdges(end)]);
plotPrefs;


%% strfs

for i = 1:length(mdls2use)
    % filter plot
    s1(i) = subplot(4,ncols,2+[0 ncols]+i);
    
    
    strf = mean(res.(mdls2use{i}).filter,3);
    clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
    cmap = flipud(cbrewer2('PiYG'));
    imagesc(res.ops.t,res.ops.f/1000,strf);
    h = colorbar; colormap(s1(i),cmap);
    caxis([clims(1) clims(end)]);
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    title(mdls2use{i});
    plotPrefs; axis square;
    
end



%% nonlinearities

for i = 1:length(mdls2use)
    subplot(4,ncols,2+[ncols*2 ncols*3]+i); hold on;

    % static model
    lbl = sprintf('%s_static',mdls2use{i});
    if isfield(res,lbl)
        tmp = vertcat(res.(lbl).xy{:});
        scatter(tmp(:,1),tmp(:,2),5,...
                'MarkerFaceColor',[.7 .7 .7],...
                'MarkerEdgeColor','none',...
                'MarkerFaceAlpha',.5)
        x = sort(tmp(:,1));
        for k = 1:size(res.(lbl).ahat,1)
            plot(x,res.model(res.(lbl).ahat(k,:),x),...
                 'Color',[.7 .7 .7],...
                 'LineWidth',1);  
        end
    end
    
    tmp = vertcat(res.(mdls2use{i}).xy{:});
    xr(i,:) = [min(tmp(:,1)) max(tmp(:,1))];
    yr(i,:) = [min(tmp(:,2)) max(tmp(:,2))];
    
    % plot nonlinearity scatter and fits
    for j = 1:2
        tmp = vertcat(res.(mdls2use{i}).xy{:,j});
        s = scatter(tmp(:,1),tmp(:,2),5);
        s.MarkerFaceColor = cc(j,:);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceAlpha = .5;
        x = sort(tmp(:,1));
        for k = 1:size(res.(mdls2use{i}).ahat,1)
            plot(x,res.model(res.(mdls2use{i}).ahat(k,j,:),x),...
                 'Color',cc(j,:),...
                 'LineWidth',1);  
        end
    end
    xlabel('Filter Response'); ylabel('FR (Hz)');
    plotPrefs; axis square;
end

xl = round(max(abs(xr(:))))+1;
yl = round(max(abs(yr(:))))+1;

for i = 1:length(mdls2use)
    s2 = subplot(4,ncols,2+[ncols*2 ncols*3]+i);
    xlim([-xl xl]); ylim([0 yl]);
    s2.Position = [s2.Position(1:2) s1(i).Position(3:4)];
end