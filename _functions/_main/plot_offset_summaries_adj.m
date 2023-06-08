function plot_offset_summaries_adj(r,ind,colors,include)

nrows = 4;
ncols = 2;
x = combindstats(r.offs,ind);
yl = [.5 1];

%% population plots
% plot behavioral pc
subplot(nrows,ncols,1); hold on;
for i = 1:2
    pc = PDtoPC(r.beh_rate_adj(r.contrastI==(i-1) & include,2,:),...
                r.beh_rate_adj(r.contrastI==(i-1) & include,1,:));
    [ym,ye] = combindstats(pc,ind);
    eh(i) = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl); 
eh(end+1) = plot([0 0],ylim,'k'); plotPrefs;
legend(eh,'H->L','L->H','switch','location','se');
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{behavior}'); axis tight; xlim([-.05 1.05]);
title('Behavior');

% plot population auc
subplot(nrows,ncols,3); hold on;
for i = 1:2
    [ym,ye] = combindstats(...
        r.auc(r.contrastI==(i-1) & include,1,:),ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl);
plot([0 0],ylim,'k'); plotPrefs;
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{population auc}'); axis tight; xlim([-.05 1.05]);
title('Population Metrics');


% plot population criterion pc
subplot(nrows,ncols,5); hold on;
for i = 1:2
    pc = PDtoPC(r.critp_adj(r.contrastI==(i-1) & include,2,:),...
                r.critp_adj(r.contrastI==(i-1) & include,1,:));
    [ym,ye] = combindstats(pc,ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl); 
plot([0 0],ylim,'k'); plotPrefs;
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{population criterion}'); axis tight; xlim([-.05 1.05]);

% plot population svm
subplot(nrows,ncols,7); hold on;
for i = 1:2
    pc = PDtoPC(r.svm_rate_adj(r.contrastI==(i-1) & include,2,:),...
                r.svm_rate_adj(r.contrastI==(i-1) & include,1,:));
    [ym,ye] = combindstats(pc,ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl); 
plot([0 0],ylim,'k'); plotPrefs;
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{svm}'); axis tight; xlim([-.05 1.05]);
xlabel('Target Offset (ms, rel. switch)');


%% single cell plots

% mean fr
lw = [.5 1];
subplot(nrows,ncols,2); hold on;
for i = 1:2
    mkfill = {'none',colors{i}{end}};
    for j = 1:2
        [ym,ye] = combindstats(...
            r.mean_fr(r.contrastI==(i-1) & include,j,:),ind);
        eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                       'Marker','o','MarkerFaceColor',mkfill{j},'LineWidth',lw(j));
    end
end
plot([0 0],ylim,'k'); plotPrefs; xlim([-.05 1.05]);
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('FR_{sig neurons}'); title('Single Cell Metrics');

% mean auc
subplot(nrows,ncols,4); hold on;
for i = 1:2
    [ym,ye] = combindstats(...
        r.mean_auc(r.contrastI==(i-1) & include,1,:),ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl); 
plot([0 0],ylim,'k'); plotPrefs; xlim([-.05 1.05]);
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{auc sig neurons}');

% mean criterion percent correct
subplot(nrows,ncols,6); hold on;
for i = 1:2
    pc = PDtoPC(r.mean_critp_adj(r.contrastI==(i-1) & include,2,:),...
                r.mean_critp_adj(r.contrastI==(i-1) & include,1,:));
    [ym,ye] = combindstats(pc,ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
%ylim(yl); 
plot([0 0],ylim,'k'); plotPrefs; xlim([-.05 1.05]);
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('PC_{criterion sig neurons}');

% prop sig neurons
subplot(nrows,ncols,8); hold on;
for i = 1:2
    [ym,ye] = combindstats(...
        r.prop_sig_neurons(r.contrastI==(i-1) & include,1,:),ind);
    eh = errorBars(x,ym,colors{i}{end},[],[],ye,[],...
                   'Marker','.','MarkerSize',20,'LineWidth',1);
end
plot([0 0],ylim,'k'); plotPrefs; axis tight; xlim([-.05 1.05]);
set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
ylabel('Proportion Significant Neurons');
xlabel('Target Offset (ms, rel. switch)');

