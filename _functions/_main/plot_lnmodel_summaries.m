function stat = plot_lnmodel_summaries(res_ln,res_psych,r_psych,ops)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model comparison for all neurons
% noise ratio cutoff
nrcut = 100;
NR = squeeze(median(cat(3,res_ln.noiseRatio),1,'omitnan'))';
I = all(NR < nrcut,2) & ~any(NR < 0,2);

% exclude all of the neurons with bad NR
rall = res_ln(I);

% contrast index [1 2] if low to high, [2 1] if high to low
nr = squeeze(median(cat(3,rall.noiseRatio),1,'omitnan'))';
cI = zeros(size(nr));
contrast = {res_ln.cond};
clh = contains(contrast(I),'lohi');
chl = ~clh;
cI(clh,:) = repmat([1 2],sum(clh),1);
cI(chl,:) = repmat([2 1],sum(chl),1);

% plots
nrows = 7;
ncols = 3;

% model comparisons
cmp = {'sta','NRC';'sta','glm';'NRC','glm'};

% results for each model
sz = [500 1100];
f1 = figure(1); set(f1,'Position',[0 0 sz]); clf;
statcnt = 0;
for i = 1:length(ops.models)

    
    %% compare across models
    % static models
    subplot(nrows,ncols,0+(i)); hold on;
    x = [rall.(sprintf('%s_static_allr',cmp{i,1}))];
    y = [rall.(sprintf('%s_static_allr',cmp{i,2}))];
    scatter(x,y,20,'k.')
    plot(median(x),median(y),'rx');
    lims = [min([x(:);y(:)]), max([x(:);y(:)])];
    plot(ylim,ylim,'r'); axis tight; plotPrefs;
    xlabel(sprintf('r (static %s)',cmp{i,1}));
    ylabel(sprintf('r (static %s)',cmp{i,2}));
    plotPrefs;
    
    % stats
    statcnt = statcnt + 1;
    [pv,~,stats] = signrank(x,y);
    stat(statcnt).type = sprintf('%s vs %s static all corr',cmp{i,1},cmp{i,2});
    stat(statcnt).test = 'signrank';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).median = [median(x) median(y)];
    stat(statcnt).iqr = [iqr(x) iqr(y)];
    stat(statcnt).n = [length(x) length(y)];
    stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
    stat(statcnt).effectSize = stats.zval / sqrt(length(x));
    
    % GC models
    subplot(nrows,ncols,ncols+(i)); hold on;
    x = [rall.(sprintf('%s_allr',cmp{i,1}))];
    y = [rall.(sprintf('%s_allr',cmp{i,2}))];
    scatter(x,y,20,'k.')
    plot(median(x),median(y),'rx');
    lims = [min([x(:);y(:)]), max([x(:);y(:)])];
    plot(ylim,ylim,'r'); axis tight; plotPrefs;
    xlabel(sprintf('r (GC %s)',cmp{i,1}));
    ylabel(sprintf('r (GC %s)',cmp{i,2}));
    plotPrefs;
    
    % stats
    statcnt = statcnt + 1;
    [pv,~,stats] = signrank(x,y);
    stat(statcnt).type = sprintf('%s vs %s GC all corr',cmp{i,1},cmp{i,2});
    stat(statcnt).test = 'signrank';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).median = [median(x) median(y)];
    stat(statcnt).iqr = [iqr(x) iqr(y)];
    stat(statcnt).n = [length(x) length(y)];
    stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
    stat(statcnt).effectSize = stats.zval / sqrt(length(x));
    
    
    %% compare within models
    % plot static versus gain control
    subplot(nrows,ncols,ncols*2+(i)); hold on
    x = [rall.(sprintf('%s_static_allr',ops.models{i}))];
    y = [rall.(sprintf('%s_allr',ops.models{i}))];
    scatter(x,y,20,'k.')
    plot(median(x),median(y),'rx');
    lims = [min([x(:);y(:)]), max([x(:);y(:)])];
    plot(ylim,ylim,'r'); axis tight; plotPrefs;
    xlabel(sprintf('r (static %s)',ops.models{i}));
    ylabel(sprintf('r (GC %s)',ops.models{i}));
    plotPrefs;
    
    % stats
    statcnt = statcnt + 1;
    [pv,~,stats] = signrank(x,y);
    stat(statcnt).type = sprintf('%s static vs GC all corr',ops.models{i});
    stat(statcnt).test = 'signrank';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).median = [median(x) median(y)];
    stat(statcnt).iqr = [iqr(x) iqr(y)];
    stat(statcnt).n = [length(x) length(y)];
    stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
    stat(statcnt).effectSize = stats.zval / sqrt(length(x));
    
    
    %% gain control
    subplot(nrows,ncols,ncols*3+i); hold on;
    p = cat(3,rall.(sprintf('%s_ahat',ops.models{i})));
    g = squeeze(p(:,3,:))'; glo = g(cI==1); ghi = g(cI==2);
    edges = linspace(-.1,.6,30);
    histogram(glo,edges,'facecolor','b');
    histogram(ghi,edges,'facecolor','r');
    ylim([0 1500]);
    plot(median([glo glo]),ylim,'b--','linewidth',1);
    plot(median([ghi ghi]),ylim,'r--','linewidth',1);
    xlabel(sprintf('Gain (%s)',ops.models{i})); 
    ylabel('Cell Count'); plotPrefs;
    
    % stats
    statcnt = statcnt + 1;
    [pv,~,stats] = signrank(glo,ghi);
    stat(statcnt).type = sprintf('%s lo vs high gain',ops.models{i});
    stat(statcnt).test = 'signrank';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).median = [median(glo) median(ghi)];
    stat(statcnt).iqr = [iqr(glo) iqr(ghi)];
    stat(statcnt).n = [length(glo) length(ghi)];
    stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
    stat(statcnt).effectSize = stats.zval / sqrt(length(x));
    
    
    %% fit in noise vs target periods
    subplot(nrows,ncols,ncols*4+i); hold on;
    x = [1 2 4 5];
    d = cat(1,rall.(sprintf('%s_meanr',ops.models{i})));
    dd{1} = d(cI(:,1)==1,1); % low contrast adapt
    dd{2} = d(cI(:,2)==1,2); % low contrast target
    dd{3} = d(cI(:,1)==2,1); % high contrast adapt
    dd{4} = d(cI(:,2)==2,2); % high contrast target
    distributionPlot(dd([1,3]),'widthDiv',[2 1],'histOri','left', ...
                     'color',{[.7 .7 1],[1 .7 .7]},'showMM',0);
    distributionPlot(dd([2,4]),'widthDiv',[2 2],'histOri','right',...
                     'color',{'b','r'},'showMM',0)
    ylabel(sprintf('r (%s)',ops.models{i})); plotPrefs;
    title('Light = adapt, dark = target');
    
    % stats (unbalanced n-way anova)
    statcnt = statcnt + 1;
    data = d(:);
    cont = cI(:);
    period = [ones(length(d),1); ones(length(d),1)*2];
    [pv,tbl,stats] = anovan(data,{cont,period},'display','off',...
                            'varnames',{'contrast','trial period'},...
                            'model',3);
    stat(statcnt).type = sprintf('contrast x trial period anova on %s model corr',...
                                 ops.models{i});
    stat(statcnt).test = 'anovan';
    stat(statcnt).p = pv;
    stat(statcnt).stats = tbl;
    stat(statcnt).n = length(d);
    stat(statcnt).effectSizeMethod = 'eta^2';
    stat(statcnt).effectSize = [tbl{2:4,2}] ./ tbl{6,2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1 2],'display','off');
    stat(statcnt).multcomp(1).comps = c;
    stat(statcnt).multcomp(1).means = m;
    stat(statcnt).multcomp(1).gnames = gnames;
    
    
    
    %% compare gain in adapt vs target periods across contrasts
    subplot(nrows,ncols,ncols*5+i); hold on;
    x = [1 2 4 5];
    [g,ind] = rmoutliers(g);
    gg{1} = g(cI(~ind,1)==1,1); % low contrast adapt
    gg{2} = g(cI(~ind,2)==1,2); % low contrast target
    gg{3} = g(cI(~ind,1)==2,1); % high contrast adapt
    gg{4} = g(cI(~ind,2)==2,2); % high contrast target
    distributionPlot(gg([1,3]),'widthDiv',[2 1],'histOri','left', ...
                     'color',{[.7 .7 1],[1 .7 .7]},'showMM',0);
    distributionPlot(gg([2,4]),'widthDiv',[2 2],'histOri','right',...
                     'color',{'b','r'},'showMM',0)
    ylim([-.1,.3]);
    ylabel(sprintf('Gain (%s)',ops.models{i})); plotPrefs;
    title('Light = adapt, dark = target');
        
    % stats (unbalanced n-way anova)
    statcnt = statcnt + 1;
    data = g(:);
    cont = cI(~ind,:);
    cont = cont(:);
    period = [ones(length(g),1); ones(length(g),1)*2];
    [pv,tbl,stats] = anovan(data,{cont,period},'display','off',...
                            'varnames',{'contrast','trial period'},...
                            'model',3);
    stat(statcnt).type = sprintf('contrast x trial period anova on %s gain',...
                                 ops.models{i});
    stat(statcnt).test = 'anovan';       
    stat(statcnt).p = pv;
    stat(statcnt).stats = tbl;
    stat(statcnt).n = length(g);
    stat(statcnt).effectSizeMethod = 'eta^2';
    stat(statcnt).effectSize = [tbl{2:4,2}] ./ tbl{6,2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1 2],'display','off');
    stat(statcnt).multcomp(1).comps = c;
    stat(statcnt).multcomp(1).means = m;
    stat(statcnt).multcomp(1).gnames = gnames;
    

end


%% noise ratio
subplot(nrows,ncols,ncols*6+1); hold on;
x = nr(cI==1);
y = nr(cI==2);
scatter(x,y,20,'k.')
plot(median(x),median(y),'rx');
lims = [min([x(:);y(:)]), max([x(:);y(:)])];
plot(ylim,ylim,'r'); axis tight; plotPrefs;
xlabel('NR (low contrast)');
ylabel('NR (high contrast)');
plotPrefs;

% stats
statcnt = statcnt + 1;
[pv,~,stats] = signrank(x,y);
stat(statcnt).type = 'low vs high contrast NR';
stat(statcnt).test = 'signrank';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = [median(x) median(y)];
stat(statcnt).iqr = [iqr(x) iqr(y)];
stat(statcnt).n = [length(x) length(y)];
stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
stat(statcnt).effectSize = stats.zval / sqrt(length(x));


saveFigPDF(f1,sz,'./_plots/_lnmodel_comparison.pdf');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship between gain and psychometric performance

%% setup

% first find the shared cells between ln model fits and
% psychometric sessions
pcells = {res_psych.single_cell.cellID};
gcells = {res_ln.cellID};
[C,ia,ib] = intersect(gcells,pcells);
rpsy = res_psych.single_cell(ib);
rln = res_ln(ia);

% noise ratio
clear NR;
NR = squeeze(median(cat(3,rln.noiseRatio),1,'omitnan'))';

% get gain and outliers
p = cat(3,rln.(sprintf('%s_ahat',ops.models{2})));
g = squeeze(p(:,3,:))';
[~,g_outliers] = rmoutliers(g);

% cell groups
target_responsive = sum(cat(1,rpsy.sig),2) > 1;
noise_responsive = all(NR < nrcut,2) & ~any(NR < 0,2);

% false alarm rate cutoff
facut = .3;

% titles
titles{1} = 'N responsive';
%titles{2} = 'T responsive';


% set up each group
index{1} = noise_responsive & ~g_outliers;
%index{2} = target_responsive & ~g_outliers;

sz2 = [350 1260];

% plot stuff
for i = 1:length(titles)
    fh = figure(i+23); clf; set(fh,'Position',[0 0 sz2]);
    title(sprintf('%s (n = %d neurons)',titles{i},sum(index{i})));
    stat = plot_gain_psych_stats(rln,r_psych,g,index{i},...
                                 stat,statcnt,facut);
    
    saveFigPDF(fh,sz2,'./_plots/_lnmodel_psych.pdf');
    
end