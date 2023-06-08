function plot_noise_clusters(res,ops)

% plot options
off = .4;
k = ops.kpcs;
lw = linspace(1.5,.25,k);
clim = .7;
nrows = 4; ncol = 4;

% all data sorted by contrast
[~,sortI] = sort(res.contrastI);
subplot(nrows,ncol,[1 5]);
colormap('default');
imagesc(ops.time,1:length(sortI),normalize(res.mPSTH(sortI,:),2),[-3 3]); colorbar;
title('z-scored FR'); plotPrefs;
ylabel('Neuron'); xlabel('Time (s)');

% high -> low
subplot(nrows,ncol,2); hold on;
plot(res.pca_explained{1}(1:10),'k');
ylabel('Variance Explained'); plotPrefs; title('Variance Explained');
cmap = [linspace(0,clim,k)' linspace(0,clim,k)' ones(k,1) ];
for i = k:-1:1
    subplot(nrows,ncol,2); hold on
    scatter(i,res.pca_explained{1}(i,:),30,cmap(i,:),'o','filled');
    
    subplot(nrows,ncol,3:4); hold on;
    plot(ops.time(ops.time_index),res.pca_coeffs{1}(:,i)+off*i,'color',cmap(i,:),'linewidth',lw(i));
end
yt = off*[1:k];
set(gca,'ytick',yt);
set(gca,'yticklabels',num2str([1:k]'));
ylabel('PC'); plotPrefs; title('PC Coefficients'); axis tight


% low -> high
subplot(nrows,ncol,6); hold on
plot(res.pca_explained{2}(1:10),'k');
ylabel('Variance Explained'); xlabel('PC'); plotPrefs; 
cmap = [ones(k,1) linspace(0,clim,k)' linspace(0,clim,k)'];
for i = k:-1:1
    subplot(nrows,ncol,6); hold on
    scatter(i,res.pca_explained{2}(i,:),30,cmap(i,:),'o','filled');
    
    subplot(nrows,ncol,7:8); hold on;
    plot(ops.time(ops.time_index),res.pca_coeffs{2}(:,i)+off*i,'color',cmap(i,:),'linewidth',lw(i));
end
yt = off*[1:k];
set(gca,'ytick',yt);
set(gca,'yticklabels',num2str([1:k]'));
ylabel('PC');
plotPrefs; axis tight;


off = 3;
for j = 1:2     

    thispsth = res.mPSTH(res.contrastI == (j-1),ops.time_index);
    
    
    % unique labels
    uL = unique(res.tsne_label{j});
    kk = length(uL);
    
    if j == 1
        cmap = [linspace(0,clim,kk)' linspace(0,clim,kk)' ones(kk,1)];
    else
        cmap = [ones(kk,1) linspace(0,clim,kk)' linspace(0,clim,kk)'];
    end
    
    
    subplot(nrows,ncol,[9 10]+4*(j-1)); hold on;
    set(gca,'ydir','normal')
    colormap('gray')
    imagesc(res.tsne_bins{j}{1},res.tsne_bins{j}{2},res.tsne_pd{j}'.*res.tsne_wslabels{j}'); colorbar;
    scatter(res.tsne_embedding{j}(:,1),res.tsne_embedding{j}(:,2),50,'k.');
    cvec = cmap(res.tsne_label{j},:);
    scatter(res.tsne_embedding{j}(:,1),res.tsne_embedding{j}(:,2),20,cvec,'.');
    axis tight;
    plotPrefs;
    if j == 1; title('tSNE Map'); end

    
    subplot(nrows,ncol,[11 12]+4*(j-1)); hold on;
    for i = 1:length(uL)
        plot(ops.time(ops.time_index),mean(normalize(thispsth(res.tsne_label{j}==i,:),2),1)+off*i,...
             'color',cmap(i,:),'linewidth',1);
    end
    yt = off*[1:kk];
    set(gca,'ytick',yt);
    set(gca,'yticklabels',num2str([1:kk]'));
    if j == 1; title('tSNE Clusters'); end
    ylabel('Cluster Number');
    plotPrefs;
    
end
xlabel('Time (s. rel. trial onset)');