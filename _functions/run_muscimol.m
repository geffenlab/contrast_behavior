function stat = run_muscimol

%addpath(genpath('~/gits/gain-gonogo/'));
%addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('./_functions/'));


%% directory stuff
dataDir = '~/gits/gain-gonogo/_data';
mouseList = {'CA123','CA124','CA125','CA126'};

% note labels
conds = {'muscimol','saline'};

clear stat;

if ~exist('./_data/_res_behavior_muscimol.mat','file');
    clear S;
    cnt = 1;
    for m = 1:length(mouseList)
        
        % mouse directory
        mdir = fullfile(dataDir,mouseList{m});
        
        % file list
        [fileList fileInd p] = indexDataFiles(mdir);
        
        % get file notes
        clear note label;
        for i = 1:length(p)
            if isfield(p{i},'note') & ~isempty(p{i}.note)
                note{i} = p{i}.note;
            else
                note{i} = ' ';
            end
            label{i} = p{i}.stimLabel;
        end
        
        fprintf('Mouse %d/%d ',m,length(mouseList)); tic;
        
        
        % for each condition (muscimol and saline control)
        for i = 1:length(conds)
            
            I = find(contains(note,conds{i})' & ...
                     contains(label,'testing')' & ...
                     ~contains(note,'discard')');
            
            for j = 1:length(I)
                
                % load data
                dat = load(fileList{I(j)});
                
                % compute performance, averaging over time
                [rate,dp,pc,nresp,ntrials] = gonogoPerformance(dat.resp,dat.tt(:,1:2),dat.abort,1,[],2);
                
                % fit percent correct with logistic function
                if size(dat.params.targetDBShift,2) == 6
                    x = dat.params.targetDBShift;
                else
                    x = dat.params.attenuation;
                end
                y = pc';
                if any(isnan(y))
                    keyboard
                end
                [prms,mdl,thresh,sense] = fitLogGrid(x,y,[],ntrials');
                
                % fit rate
                y = rate(2:end)';
                [rate_prms,mdl,rate_thresh] = fitLogGrid(x,y,[],ntrials');
                
                % performance at threshold
                rate_at_thresh = mdl(rate_prms,thresh);
                
                % probability of licking over all trials
                [fp,fn] = fileparts(fileList{I(j)});
                [t,trialType,response,RT,abort] = ...
                    parseLog(fullfile(fp,[fn '.txt']));
                edges = 1e6*(-.1:.025:5.2);
                lickPSTH = makePSTH(vertcat(t.lickTimes),vertcat(t.stimTimes),edges);
                lickPSTH = lickPSTH * mean(diff(edges));
                
                
                % save to struct
                S(cnt).mouse = mouseList{m};
                S(cnt).session = dat.params.IDsess;
                S(cnt).contrast = dat.params.sd(2) > dat.params.sd(1);
                S(cnt).condition = conds{i};
                S(cnt).notes = dat.params.note;
                S(cnt).snr = x;
                S(cnt).rate = rate';
                S(cnt).dp = dp';
                S(cnt).pc = pc';
                S(cnt).prms = prms';
                S(cnt).thresh = thresh;
                S(cnt).sense = sense;
                S(cnt).mxslope = max(diff(pc')./diff(x));
                S(cnt).rate_prms = rate_prms';
                S(cnt).rate_thresh = rate_thresh;
                S(cnt).rate_at_thresh = rate_at_thresh;
                S(cnt).rate_mxslope = max(diff(rate(2:end)')./diff(x));
                S(cnt).plick = sum(lickPSTH) / size(lickPSTH,1);
                S(cnt).lickcnt = sum(lickPSTH);
                S(cnt).label = label{I(j)};
                
                cnt = cnt + 1;
                
            end
            
        end
        
        toc;
        
    end

    T = struct2table(S);
    T.fa = T.rate(:,1);

    lickTime = edges / 1e6;
    lickTime = lickTime(1:end-1) - mean(diff(lickTime))/2;

    T = T(~contains(T.label,'NoNoise'),:);
    save('./_data/_res_behavior_muscimol.mat');

else
    load('./_data/_res_behavior_muscimol.mat');
end


%% by mouse
groupVars = {'mouse','condition','contrast'};
statVars = {'snr','pc','rate','fa','thresh','mxslope','maxrate','plick'};
T.maxrate = T.rate(:,end);
m = grpstats(T,groupVars,'mean','datavars',statVars);

% mouse average stats (based on average curve, not session means)
for i = 1:size(m,1)
    
    I = contains(T.mouse,m.mouse{i})' & ...
        contains(T.condition,m.condition{i})' & ...
        T.contrast' == m.contrast(i);
    
    % fit pc mean
    x = T.snr(I,:); y = T.pc(I,:);
    [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
    slope = max(diff(mean(y,1))./diff(mean(x,1)));
    
    m.thresh(i) = thresh;
    m.slope(i) = slope;
    m.params(i,:) = prms';
    
    % fit rate mean
    y = T.rate(I,2:end);
    [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
    slope = max(diff(mean(y,1))./diff(mean(x,1)));
    
    m.thresh_rate(i) = thresh;
    m.slope_rate(i) = slope;
    m.params_rate(i,:) = prms';

end


% plot vars
colors = {'b','r'};
cc = [0 0 1; 1 0 0];
cclite = [.7 .7 1; 1 .7 .7];
alphas = [.2 .6];
ms = 6;
ls = {'--','-'};
xf = linspace(0,25,100);


%% individual mice
f0 = figure(10); clf;
% plot mouse averages and fits
for i = 1:length(mouseList)
    
    % for each contrast
    for j = 1:2
        
        I = contains(m.mouse,mouseList{i}) & ...
            m.contrast == (j-1);
        mfill = {'none',cclite(j,:)};

        
        if sum(I) > 1
            plotI = sub2ind([2 4],j,i);
            subplot(4,2,plotI); hold on
            
            % plot muscimol and saline conditions
            for k = 1:2
                
                ind = contains(T.mouse,mouseList{i}) & ...
                      T.contrast == (j-1) & ...
                      contains(T.condition,conds{k});
                mI = contains(m.mouse,mouseList{i}) & ...
                     m.contrast == (j-1) & ...
                     contains(m.condition,conds{k});
                
                % plot mean data
                x = mean(T.snr(ind,:),1);
                ym = mean(T.rate(ind,2:end),1);
                ye = sem(T.rate(ind,2:end),1);
                errorbar(x,ym,ye,'o','color',cclite(j,:),'markerfacecolor',...
                         mfill{k},'markersize',ms);
                errorbar(-5,mean(T.rate(ind,1),1),sem(T.rate(ind,1),1),'o',...
                         'color',cclite(j,:),'markerfacecolor',mfill{k},...
                         'markersize',ms);
                
                % plot rate fits
                plot(xf,mdl(m.params_rate(mI,:),xf),'color',colors{j},...
                     'linestyle',ls{k},'linewidth',1);
                title(mouseList{i}); plotPrefs;
                
                if (i+j+k) == 8
                    set(gca,'Xtick',[-5 x],'xticklabel',['FA'; num2str(x')]);
                    ylabel('p(Respond)'); 
                    xlabel('Target Volume (dB SNR)');
                end
            end
        end
    end
end

saveFigPDF(f0,[400 600],'./_plots/_muscimol_psych_individual.pdf');













%% stats across sessions
sz1 = [400 1200];
f1 = figure(1); clf; set(f1,'Position', [0 0 sz1]);
clear p;

% stats
groupVars = {'session','mouse','condition','contrast'};
T.maxrate = T.rate(:,end);
T.cumLickP = cumsum(T.plick,2) ./ sum(T.plick,2);
statVars = {'snr','pc','rate','fa',...
            'thresh','mxslope','maxrate','rate_thresh','rate_mxslope',...
            'plick','cumLickP'};
m = grpstats(T,groupVars,'mean','datavars',statVars);
tableInfo(T);

% plot params
% plot vars
colors = {'b','r'};
cc = [0 0 1; 1 0 0];
cclite = [.7 .7 1; 1 .7 .7];
alphas = [.2 .6];
ms = 6;
ls = {'--','-'};
xf = linspace(0,25,100);
nrows = 6;
ncols = 2;
contrastLabel = {'low contrast','high contrast'};

% initialize structure for statistics
stat = struct(); statcnt = 0;

% for each contrast
for j = 1:2
    mfill = {'none',cc(j,:)};
    subplot(nrows,ncols,j); hold on;
    for i = 1:length(conds)
        
        I = contains(m.condition,conds{i}) & m.contrast == (j-1);
        x = m.mean_snr(I,:); xf = linspace(min(x(:)),max(x(:)),100);
        y = m.mean_rate(I,2:end);
        [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
        
        plot(mean(x,1)',y','color',cclite(j,:),'linestyle',ls{i}),
        errorbar(mean(x,1),mean(y,1),sem(y,1),'o','color',cc(j,:),...
                 'markerfacecolor',mfill{i},'markersize',ms);
        p(i) =  plot(xf,mdl(prms,xf),'color',colors{j},'linestyle',ls{i},...
                     'linewidth',1);
        errorbar(min(x(:)) - 5,mean(m.mean_rate(I,1),1),sem(m.mean_rate(I,1),1),...
                 'o','color',cc(j,:),'markerfacecolor',...
                 mfill{i},'markersize',ms);
        
    end
    legend(p,conds,'location','nw');
    set(gca,'xtick',[-5 mean(x)])
    set(gca,'xticklabels',num2str([-inf mean(x)]'))
    xlabel('Target Volume (dB SNR)'); xlim([-10 30]);
    ylabel('P(respond)'); ylim([0 1]);
    plotPrefs;

end

stats2run = {'mean_rate_thresh','mean_rate_mxslope','mean_maxrate','mean_fa'};
clear s;
for k = 1:length(stats2run)
    
    hold on;
    cnt = 1;
    
    yall = [];
    musc = [];
    cont = [];
    
    % for each contrast
    for i = 1:2
        
        subplot(nrows,ncols,3+(k-1)); hold on;
        xs = []; ys = [];
        
        % for each condition
        for j = 1:2
            
            I = contains(m.condition,conds{j}) & m.contrast == (i-1);            
            x = cnt; y = mean(m.(stats2run{k})(I));
            s = plotSpread(m.(stats2run{k})(I),'xValues',cnt);
            set(s{1}(1),'Marker','o','MarkerFaceColor','w',...
                        'MarkerEdgeColor','k','MarkerSize',ms-1);
            ys{j} = m.(stats2run{k})(I);
            xs = [xs cnt];            

            b(cnt) = bar(x,y);
            b(cnt).FaceColor = cc(i,:); b(cnt).FaceAlpha = alphas(j);
            b(cnt).LineStyle = ls{j}; b(cnt).EdgeColor = cc(i,:);
            b(cnt).LineWidth = 1;
            
            cnt = cnt + 1;
            if i == 1 & j == 2
                cnt = cnt + 1;
            end
            
            % aggregates
            yall = [yall; ys{j}];
            musc = [musc; repmat(j,length(ys{j}),1)];
            cont = [cont; repmat(i,length(ys{j}),1)];

        end
        
        
        % plot lines and pval d
        [pv,~,stats] = ranksum(ys{1},ys{2});
        plotPvalSmart(pv,xs,max(cellfun(@max,ys)));
        
        % stats
        statcnt = statcnt + 1;
        stat(statcnt).type = sprintf('%s %s muscimol vs saline noise',stats2run{k},contrastLabel{i});
        stat(statcnt).test = 'ranksum';
        stat(statcnt).p = pv;
        stat(statcnt).stats = stats;
        stat(statcnt).median = cellfun(@median,ys);
        stat(statcnt).iqr = cellfun(@iqr,ys);
        stat(statcnt).n = cellfun(@numel,ys);
        stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
        stat(statcnt).effectSize = stats.zval / sqrt(sum(cellfun(@numel,ys)));        
        
    end
    
    % n way anova on this statistic
    statcnt = statcnt + 1;
    [pv,tbl,stats] = anovan(yall,{musc,cont},'model','interaction',...
                            'varnames',{'muscimol/saline','contrast'},...
                            'display','off');
    [comp,means,~,gnames] = multcompare(stats,'dimension',[1 2],'display','off');
    stat(statcnt).type = sprintf('%s muscimolxcontrast anova',stats2run{k});
    stat(statcnt).test = 'anovan';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).table = tbl;
    stat(statcnt).mean = means;
    stat(statcnt).effectSizeMethod = 'partial eta^2';
    for i = 1:3
        stat(statcnt).effectSize(i) = tbl{i+1,2} / (tbl{i+1,2} + tbl{end-1,2});
    end
    stat(statcnt).multcomp = comp;
    stat(statcnt).multcompGrps = gnames;
    
    set(gca,'xtick',[1 2 4 5]);
    set(gca,'xticklabel',{'Muscimol','Saline','Muscimol','Saline'});
    ylabel(stats2run{k},'Interpreter','none'); plotPrefs;
    xtickangle(45);
    if k > 2
        ylim([0 1]);
    end
end

% 2x2x7 ANOVA: intervention (musc/saline) x task (hi/lo) x volume
statcnt = statcnt + 1;
data = m.mean_rate;
intervention = repmat(contains(m.condition,'muscimol'),1,7);
task = repmat(m.contrast,1,7);
volume = repmat([0 1 2 3 4 5 6],size(data,1),1);
[pv,tbl,stats] = anovan(data(:),{intervention(:),task(:),volume(:)},...
                        'varnames',{'intervention','task','volume'},...
                        'display','off','model','interaction');
stat(statcnt).type = '3 way ANOVA: musc/saline x hi/lo contrast x volume ~ hit rate';
stat(statcnt).test = 'anovan';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).table = tbl;
stat(statcnt).effectSizeMethod = 'partial eta^2';
for i = 1:6
    stat(statcnt).effectSize(i) = tbl{i+1,2} / (tbl{i+1,2} + tbl{end-1,2});
end
for i = 1:2
    [comp,means,~,gnames] = multcompare(stats,'dimension',i,'display','off');
    stat(statcnt).multcomp(i,:) = comp;
    stat(statcnt).multcompGrps{i} = gnames;
    stat(statcnt).mean{i} = means;
end

%%
% for each contrast plot percent correct performance (not raw)
for j = 1:2
    mfill = {'none',cc(j,:)};
    subplot(nrows,ncols,6+j); hold on;
    plot([-5 30],[.5 .5],'k');
    for i = 1:length(conds)
        
        I = contains(m.condition,conds{i}) & m.contrast == (j-1);
        x = m.mean_snr(I,:); xf = linspace(min(x(:)),max(x(:)),100);
        y = m.mean_pc(I,:);
        [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
        
        plot(mean(x,1)',y','color',cclite(j,:),'linestyle',ls{i});
        errorbar(mean(x,1),mean(y,1),sem(y,1),'o','color',cc(j,:),...
                 'markerfacecolor',mfill{i},'markersize',ms);
        p(i) = plot(xf,mdl(prms,xf),'color',colors{j},'linestyle',ls{i},...
                    'linewidth',1);
        
    end
    set(gca,'xtick',[mean(x)])
    set(gca,'xticklabels',num2str([mean(x)]'))
    xlabel('Target Volume (dB SNR)'); xlim([-5 30]);
    ylabel('PC_{behavior}'); ylim([.4 1]);
    plotPrefs;

end


%%
% plot lick probability with and without muscimol
s9 = subplot(nrows,ncols,9); hold on;
mcolors = {'g','k'};
pp(1) = patchErrorBars(lickTime,m.mean_plick(contains(m.condition,'muscimol'),:),...
                       mcolors{1},[],2);
pp(2) = patchErrorBars(lickTime,m.mean_plick(contains(m.condition,'saline'),:),...
                       mcolors{2},[],2);
xlabel('Lick Time (s, rel. trial start)'); axis tight;
ylabel('Lick Probability'); legend(pp,'muscimol','saline','location','nw');
plotPrefs;

% cumulative probability
s10 = subplot(nrows,ncols,10); hold on;
pp(1) = patchErrorBars(lickTime,m.mean_cumLickP(contains(m.condition,'muscimol'),:),...
                       mcolors{1},[],2);
pp(2) = patchErrorBars(lickTime,m.mean_cumLickP(contains(m.condition,'saline'),:),...
                       mcolors{2},[],2);
xlabel('Lick Time (s, rel. trial start)');
ylabel('Cumulative Probability'); plotPrefs;
s10.Position(1) = s10.Position(1) + .1;
s10.Position(3) = s10.Position(3) * .7;
s9.Position(3) = s9.Position(3) * 1.3;

% early time window, look at lick probability
s11 = subplot(nrows,ncols,11); hold on;
tI = lickTime > 0 & lickTime < 2;
mean_lickP_early = mean(m.mean_plick(:,tI),2);
clear xs ys x y;
for i = 1:2
    
    % plots
    I = contains(m.condition,conds{i});
    x = i;
    s = plotSpread(mean_lickP_early(I),'XValues',x);
    set(s{1}(1),'Marker','o','MarkerFaceColor','w',...
                'MarkerEdgeColor','k','MarkerSize',ms-1);
    y = mean(mean_lickP_early(I));
    b(cnt) = bar(x,y);
    b(cnt).FaceColor = mcolors{i}; b(cnt).FaceAlpha = alphas(i);
    b(cnt).LineStyle = ls{i}; b(cnt).EdgeColor = mcolors{i};
    b(cnt).LineWidth = 1;
    
    % stats
    xs(i) = i;
    ys{i} = mean_lickP_early(I);
    
    
end
set(gca,'xtick',[1 2]); set(gca,'xticklabel',{'Muscimol','Saline'});
ylabel('p(lick)_{early}'); ylim([0 .06]); plotPrefs;
s11.Position(3) = s11.Position(3) * .5;


% stats
statcnt = statcnt + 1;
[pv,~,stats] = ranksum(ys{1},ys{2});
plotPvalSmart(pv,xs,max(cellfun(@max,ys)));
stat(statcnt).type = 'plick_early muscimol vs saline';
stat(statcnt).test = 'ranksum';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = cellfun(@median,ys);
stat(statcnt).iqr = cellfun(@iqr,ys);
stat(statcnt).n = cellfun(@numel,ys);
stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
stat(statcnt).effectSize = stats.zval / sqrt(sum(cellfun(@numel,ys)));  


% late time window, look at lick probability
s12 = subplot(nrows,ncols,12); hold on;
tI = lickTime > 3 & lickTime < 5;
mean_lickP_early = mean(m.mean_plick(:,tI),2);
clear xs ys x y;
for i = 1:2
    
    % plots
    I = contains(m.condition,conds{i});
    x = i;
    s = plotSpread(mean_lickP_early(I),'XValues',x);
    set(s{1}(1),'Marker','o','MarkerFaceColor','w',...
                'MarkerEdgeColor','k','MarkerSize',ms-1);
    y = mean(mean_lickP_early(I));
    b(cnt) = bar(x,y);
    b(cnt).FaceColor = mcolors{i}; b(cnt).FaceAlpha = alphas(i);
    b(cnt).LineStyle = ls{i}; b(cnt).EdgeColor = mcolors{i};
    b(cnt).LineWidth = 1;
    
    % stats
    xs(i) = i;
    ys{i} = mean_lickP_early(I);
    
    
end
set(gca,'xtick',[1 2]); set(gca,'xticklabel',{'Muscimol','Saline'});
ylabel('p(lick)_{late}'); ylim([0 .06]); plotPrefs;
s12.Position(3) = s12.Position(3) * .5;

% stats
statcnt = statcnt + 1;
[pv,~,stats] = ranksum(ys{1},ys{2});
plotPvalSmart(pv,xs,max(cellfun(@max,ys)));
stat(statcnt).type = 'plick_early muscimol vs saline';
stat(statcnt).test = 'ranksum';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = cellfun(@median,ys);
stat(statcnt).iqr = cellfun(@iqr,ys);
stat(statcnt).n = cellfun(@numel,ys);
stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
stat(statcnt).effectSize = stats.zval / sqrt(sum(cellfun(@numel,ys)));
    
saveFigPDF(f1,sz1,'./_plots/_muscimol_session_avg.pdf');








%% actx data
sz22 = [400 800];
f22 = figure(22); clf; set(f22,'Position',[0 0 sz22]);

% conditions for each mouse
mList = {'AW142','AW146'};
wd = '~/chris-lab/projects/contrast_behavior';
datDir = [wd '/_data'];
resDir = [wd '/_data'];
condLong = {'0.25mg/mL muscimol',...
            'saline'};
condShort = {'musc25','saline'}; 

% plot params
ls = {'-','--'};
ms = 6;
mstyle = {'o','o'};
colors = {'b','r'};
snrs = mean(m.mean_snr);
xtick = [-5 snrs];
xtickl = num2str([-inf snrs]');
xl = [xtick(1) xtick(end)] + mean(diff(xtick))/2 .* [-1 1];
nrows = 4;
ncols = 2;

ax(1) = subplot(nrows,ncols,1);
ax(2) = subplot(nrows,ncols,2);
ax(3) = subplot(nrows,ncols,3);
ax(4) = subplot(nrows,ncols,4);

% load cortical data
load(fullfile(datDir,'muscimolCortexData.mat'));


% run saline, then muscimol single cell analysis
order = [2 1]; clear targMean targMeanPre targAUC
if ~exist(fullfile(datDir,'_res_cortex_muscimol.mat'))
    for i = 1:length(order)
        
        I = order(i);
        % run single cell analysis
        tic;
        [~,rr{i}] = muscimol_res(musc_ctx_data(I));
        toc;
    end
    
    save(fullfile(datDir,'_res_cortex_muscimol.mat'),'rr');
else
    load(fullfile(datDir,'_res_cortex_muscimol.mat'));
end

% plot results
for i = 1:length(order)
    
    I = order(i);  
    res = rr{i};
    
    % for each cell
    for c = 1:length(res)
        
        % pre-post FR
        targMeanPre{i}(c,:,:) = squeeze(mean(res(c).condTarget{1},1));
        tmp = squeeze(mean(res(c).condTarget{2},1));
        targMean{i}(c,:,:) = tmp;
        
        % auc
        targAUCPre{i}(c,:,:) = res(c).auc{1};
        targAUCPost{i}(c,:,:) = res(c).auc{2};

    end
    
    % plot each contrast
    clear x y;
    for j = 1:2
        
        %% spks/s
        axes(ax(j)); hold on;
        h = errorBars(snrs,squeeze(targMean{i}(:,j,2:end)),cc(j,:));
        h.Marker = mstyle{I}; h.MarkerSize = ms; h.LineStyle = ls{i};
        if I == 2
            h.MarkerFaceColor = cc(j,:);
        end
        h = errorBars(xtick(1),squeeze(targMean{i}(:,j,1)),cclite(j,:));
        h.Marker = mstyle{I}; h.MarkerSize = ms;
        if I == 2
            h.MarkerFaceColor = cclite(j,:);
        end
        axis tight; plotPrefs; xlim(xl);
        set(ax(j),'xtick',xtick,'xticklabels',xtickl);
        xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
        
        %% auc
        axes(ax(j+2)); hold on;
        if i == 1
            plot(xl,[.5 .5],'k');
        end
        h = errorBars(snrs,squeeze(targAUCPost{i}(:,j,:)),cc(j,:));
        h.Marker = mstyle{I}; h.MarkerSize = ms; h.LineStyle = ls{i};
        if I == 2
            h.MarkerFaceColor = cc(j,:);
        end
        axis tight; plotPrefs; xlim(xl);
        set(ax(j+2),'xtick',xtick,'xticklabels',xtickl);
        xlabel('Target Volume (dB SNR)'); ylabel('AUC');        
    end
    
end

%% fr
% muscimol pre post
subplot(nrows,ncols,5); hold on;
I = 1; clear h;
for j = 1:2
    h(2) = errorBars(snrs+(j-1)*.5,squeeze(targMean{order(I)}(:,j,2:end)),cclite(j,:),...
              [],[],[],'median','LineStyle','--','Marker',mstyle{order(1)});
    errorBars(xtick(1)+(j-1)*.5,squeeze(targMean{order(I)}(:,j,1)),cclite(j,:),...
              [],[],[],'median','LineStyle','--','Marker',mstyle{order(1)});
    h(1) = errorBars(snrs+(j-1)*.5,squeeze(targMeanPre{order(I)}(:,j,2:end)),cc(j,:),...
              [],[],[],'median','LineStyle','--','Marker', mstyle{order(1)});
    errorBars(xtick(1)+(j-1)*.5,squeeze(targMeanPre{order(I)}(:,j,1)),cc(j,:),...
              [],[],[],'median','LineStyle','--','Marker',mstyle{order(1)});
end
title('Muscimol Pre-Post'); xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
legend(h,'Pre','Post','location','nw');  xlim(xl);
set(gca,'xtick',xtick,'xticklabels',xtickl); plotPrefs;

% saline pre post
subplot(nrows,ncols,6); hold on;
I = 2; clear h;
for j = 1:2
    h(2) = errorBars(snrs-.5,squeeze(targMean{order(I)}(:,j,2:end)),cclite(j,:),...
              [],[],[],'median','LineStyle','-','Marker',mstyle{order(1)},...
              'MarkerFaceColor',cclite(j,:));
    errorBars(xtick(1)-.5,squeeze(targMean{order(I)}(:,j,1)),cclite(j,:),...
              [],[],[],'median','LineStyle','-','Marker',mstyle{order(1)},...
              'MarkerFaceColor',cclite(j,:));
    h(1) = errorBars(snrs+.5,squeeze(targMeanPre{order(I)}(:,j,2:end)),cc(j,:),...
              [],[],[],'median','LineStyle','-','Marker', mstyle{order(1)},...
              'MarkerFaceColor',cc(j,:));
    errorBars(xtick(1)+.5,squeeze(targMeanPre{order(I)}(:,j,1)),cc(j,:),...
              [],[],[],'median','LineStyle','-','Marker',mstyle{order(1)},...
              'MarkerFaceColor',cc(j,:));
end
title('Saline Pre-Post'); xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
legend(h,'Pre','Post','location','nw'); plotPrefs;
set(gca,'xtick',xtick,'xticklabels',xtickl);  xlim(xl);

%% auc
% muscimol pre post
subplot(nrows,ncols,7); hold on;
plot([-5 30],[.5 .5],'k-');
I = 1; clear h;
for j = 1:2
    h(2) = errorBars(snrs+(j-1)*.5,squeeze(targAUCPost{order(I)}(:,j,:)),cclite(j,:),...
              [],[],[],'median','LineStyle','--','Marker',mstyle{order(1)});
    h(1) = errorBars(snrs+(j-1)*.5,squeeze(targAUCPre{order(I)}(:,j,:)),cc(j,:),...
              [],[],[],'median','LineStyle','--','Marker', mstyle{order(1)});
end
title('Muscimol Pre-Post'); xlabel('Target Volume (dB SNR)'); ylabel('AUC');
legend(h,'Pre','Post','location','nw'); ylim([.45 1]);
set(gca,'xtick',xtick,'xticklabels',xtickl); axis tight; plotPrefs;

% saline pre post
subplot(nrows,ncols,8); hold on;
plot([-5 30],[.5 .5],'k-');
I = 2; clear h;
for j = 1:2
    h(2) = errorBars(snrs-.5,squeeze(targAUCPost{order(I)}(:,j,:)),cclite(j,:),...
              [],[],[],'median','LineStyle','-','Marker',mstyle{order(1)},...
              'MarkerFaceColor',cclite(j,:));
    h(1) = errorBars(snrs+.5,squeeze(targAUCPre{order(I)}(:,j,:)),cc(j,:),...
              [],[],[],'median','LineStyle','-','Marker', mstyle{order(1)},...
              'MarkerFaceColor',cc(j,:));
end
title('Saline Pre-Post'); xlabel('Target Volume (dB SNR)'); ylabel('AUC');
legend(h,'Pre','Post','location','nw'); ylim([.45 1]);
set(gca,'xtick',snrs); axis tight; plotPrefs;

% pre-post statistics for saline
data_post = targMean{1}; data_pre = targMeanPre{1};
cgrp = repmat([ones(1,size(data_post,3)); ones(1,size(data_post,3))*2],size(data_post,1),1);
lgrp = repmat(repmat(1:7,size(data_post,2),1),size(data_post,1),1);
pre = cellstr(repmat('pre',numel(data_pre),1));
post = cellstr(repmat('post',numel(data_post),1));

y = [data_pre(:);data_post(:)];
contrast = [cgrp(:);cgrp(:)];
volume = [lgrp(:);lgrp(:)];
pre_post = [pre;post];
[pv,table,stats,terms] = anovan(y,{pre_post,contrast,volume},'display','off',...
           'model',1,'varnames',{'pre_post','contrast','volume'});
statcnt = statcnt + 1;
stat(statcnt).type = '3way anova saline pre_post, contrast, volume'
stat(statcnt).test = 'anovan';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).table = table;
stat(statcnt).n = [size(targMeanPre{1},1) size(targMean{1},1)];
stat(statcnt).effectSizeMethod = 'eta^2';
stat(statcnt).effectSize(1) = table{2,2} / table{end,2};
stat(statcnt).effectSize(2) = table{3,2} / table{end,2};
stat(statcnt).effectSize(3) = table{4,2} / table{end,2};


% pre-post statistics for muscimol
data_post = targMean{2}; data_pre = targMeanPre{2};
cgrp = repmat([ones(1,size(data_post,3)); ones(1,size(data_post,3))*2],size(data_post,1),1);
lgrp = repmat(repmat(1:7,size(data_post,2),1),size(data_post,1),1);
pre = cellstr(repmat('pre',numel(data_pre),1));
post = cellstr(repmat('post',numel(data_post),1));

y = [data_pre(:);data_post(:)];
contrast = [cgrp(:);cgrp(:)];
volume = [lgrp(:);lgrp(:)];
pre_post = [pre;post];
[pv,table,stats,terms] = anovan(y,{pre_post,contrast,volume},'display','off',...
           'model',1,'varnames',{'pre_post','contrast','volume'});
statcnt = statcnt + 1;
stat(statcnt).type = '3way anova muscimol 25ul/g pre_post, contrast, volume'
stat(statcnt).test = 'anovan';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).table = table;
stat(statcnt).n = [size(targMeanPre{2},1) size(targMean{2},1)];
stat(statcnt).effectSize(1) = table{2,2} / table{end,2};
stat(statcnt).effectSize(2) = table{3,2} / table{end,2};
stat(statcnt).effectSize(3) = table{4,2} / table{end,2};


% anova for effect of saline, level, contrast during post
data_saline = targAUCPost{1}; 
cgrp_saline = repmat([ones(1,size(data_saline,3)); ones(1,size(data_saline,3))*2],...
                     size(data_saline,1),1);
lgrp_saline = repmat(repmat(1:6,size(data_saline,2),1),size(data_saline,1),1);
data_muscimol = targAUCPost{2};
cgrp_muscimol = repmat([ones(1,size(data_muscimol,3)); ones(1,size(data_muscimol,3))*2],...
                     size(data_muscimol,1),1);
lgrp_muscimol = repmat(repmat(1:6,size(data_muscimol,2),1),size(data_muscimol,1),1);
saline = cellstr(repmat('saline',numel(data_saline),1));
muscimol = cellstr(repmat('muscimol',numel(data_muscimol),1));

y = [data_saline(:);data_muscimol(:)];
saline_muscimol = [saline;muscimol];
contrast = [cgrp_saline(:);cgrp_muscimol(:)];
volume = [lgrp_saline(:);lgrp_muscimol(:)];
[pv,table,stats,terms] = anovan(y,{saline_muscimol,contrast,volume},'display','off',...
                                'model',1,'varnames',{'saline_muscimol','contrast','volume'});
statcnt = statcnt + 1;
stat(statcnt).type = '3way anova saline_muscimol, contrast, volume'
stat(statcnt).test = 'anovan';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).table = table;
stat(statcnt).n = [size(targMean{1},1) size(targMean{2},1)];
stat(statcnt).effectSize(1) = table{2,2} / table{end,2};
stat(statcnt).effectSize(2) = table{3,2} / table{end,2};
stat(statcnt).effectSize(3) = table{4,2} / table{end,2};

saveFigPDF(f22,sz22,'./_plots/_muscimol_actx.pdf');







%% control analysis
T = struct2table(S);
T = T(contains(T.label,'NoNoise'),:);
T.fa = T.rate(:,1);
groupVars = {'session','mouse','condition'}
T.maxrate = T.rate(:,end);
T.cumLickP = cumsum(T.plick,2) ./ sum(T.plick,2);
statVars = {'snr','pc','rate','fa',...
            'thresh','mxslope','maxrate','rate_thresh','rate_mxslope',...
            'plick','cumLickP'};
m = grpstats(T,groupVars,'mean','datavars',statVars);
snrs = mean(m.mean_snr);
linestyle = {'--','-'};
facecolor = {'none','k'};

% unique mice and conditions
uM = unique(m.mouse); uC = unique(m.condition);

%  f3 = figure(3); clf;
%  % for each mouse, plot the mean and fit in each condition
% 
%  for i = 1:length(uM)
%      
%      subplot(2,1,i); hold on;
%      
%      for j = 1:length(uC)
%      
%          I = contains(m.mouse,uM{i}) & contains(m.condition,uC{j});
%          
%          % plot data
%          errorbar(snrs,mean(m.mean_rate(I,2:end)),sem(m.mean_rate(I,2:end)),...
%                   'color','k','linestyle','none','MarkerFaceColor',facecolor{j},...
%                   'MarkerSize',ms,'Marker','o');
%          errorbar(-90,mean(m.mean_rate(I,1)),sem(m.mean_rate(I,1)),...
%                   'color','k','linestyle','none','MarkerFaceColor',facecolor{j},...
%                   'MarkerSize',ms,'Marker','o');
%          
%          % fit the data
%          x = m.mean_snr(I,:); y = m.mean_rate(I,2:end);
%          [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
%          xf = linspace(min(snrs),max(snrs),100);
%          plot(xf,mdl(prms,xf),'k','LineStyle',linestyle{j});
%          xlabel(sprintf('Target Attenuation\n(dB rel. 25dB SNR)'));
%          ylabel('Response Rate');
%          set(gca,'xtick',-90:15:0);
%          set(gca,'xticklabels',num2str([-inf; snrs']));
%          plotPrefs;
%          
%      end
%  end
%  
%  saveFigPDF(f3,[180 350],'./_plots/_muscimol_control_mouse_avg.pdf');










% get the data just for these mice across all conditions in high contrast
mice = {'CA124','CA125'};
T = struct2table(S);
T = T(contains(T.mouse,mice)&T.contrast==1,:);

% nonoise conditions
nnI = contains(T.label,'NoNoise');
T.control = nnI;

f4 = figure(4); clf; nplotRows = 3;
ncol = [1 0 0; 0 0 0];
ncol_lite = [1 .7 .7; .7 .7 .7];
ls = {'--','-'};
xl = {'Target Volume (dB SNR)','Target Attenuation (dB rel 25dB SNR)'};
% plot curves for noise and nonoise muscmimol/saline and 
for i = 1:2
    
    subplot(nplotRows,2,i); hold on;
    fcol = {'none',ncol(i,:)};
    
    for j = 1:length(uC)
        
        I = T.control == (i-1) & contains(T.condition,uC{j});
        x = T.snr(I,:); y = T.rate(I,:);
        mx = mean(x,1); my = mean(y(:,2:end),1);
        xf = linspace(min(mx),max(mx),100);
        [prms,mdl] = fitLogGrid(mx,my);
        
        plot(x',y(:,2:end)','color',ncol_lite(i,:),'linestyle',ls{j});
        errorbar(mx,mean(y(:,2:end),1),sem(y(:,2:end)),...
                 'color',ncol(i,:),'linestyle','none',...
                 'marker','o','markerfacecolor',fcol{j});
        errorbar(mx(1) - mean(diff(mx)),mean(y(:,1),1),sem(y(:,1)),...
                 'color',ncol(i,:),'linestyle','none',...
                 'marker','o','markerfacecolor',fcol{j});
        plot(xf,mdl(prms,xf),'color',ncol(i,:),'linestyle',ls{j});
        set(gca,'xtick',[mx(1) - mean(diff(mx)) mx]);
        set(gca,'xticklabels',num2str([-inf;mx']));
        xlabel(xl{i}); ylabel('Response Rate');
        plotPrefs;

    end
end


T.fa = T.rate(:,1);
T.maxrate = T.rate(:,end);
T.maxrate2 = mean(T.rate(:,end-1:end),2);
stats2run = {'rate_at_thresh','rate_mxslope','maxrate','fa'};
backLbl = {'Noise','No Noise'};
condLbl = {'Muscimol','Saline'};
for i = 1:length(stats2run)
    
    subplot(nplotRows,2,2+i); hold on;
    cnt = 1;
    
    yall = [];
    musc = [];
    task = [];
    
    % for each background condition
    for j = 1:2
        
        xs = []; ys = [];
        
        % for muscimol/saline
        for k = 1:length(uC)
            
            
            I = contains(T.condition,conds{k}) & T.control == (j-1);            
            x = cnt; y = mean(T.(stats2run{i})(I));
            s = plotSpread(T.(stats2run{i})(I),'xValues',cnt);
            set(s{1}(1),'Marker','o','MarkerFaceColor','w',...
                        'MarkerEdgeColor','k','MarkerSize',ms-1);
            ys{k} = T.(stats2run{i})(I);
            xs = [xs cnt];
            
            b(cnt) = bar(x,y);
            b(cnt).FaceColor = ncol(j,:); b(cnt).FaceAlpha = alphas(k);
            b(cnt).LineStyle = ls{k}; b(cnt).EdgeColor = ncol(j,:);
            b(cnt).LineWidth = 1;
            
            cnt = cnt + 1;
            if j == 1 & k == 2
                cnt = cnt + 1;
            end
            
            % add to aggregates
            yall = [yall; ys{k}];
            musc = [musc; repmat(k,length(ys{k}),1)];
            task = [task; repmat(j,length(ys{k}),1)];
            
            
        end
        
         % plot lines and pval d
        [pv,~,stats] = ranksum(ys{1},ys{2});
        plotPvalSmart(pv,xs,max(cellfun(@max,ys)));
        statcnt = statcnt + 1;
        stat(statcnt).type = sprintf('%s %s muscimol vs saline',stats2run{i},backLbl{j});
        stat(statcnt).test = 'ranksum';
        stat(statcnt).p = pv;
        stat(statcnt).stats = stats;
        stat(statcnt).median = cellfun(@median,ys);
        stat(statcnt).iqr = cellfun(@iqr,ys);
        stat(statcnt).n = cellfun(@numel,ys);
        stat(statcnt).effectSizeMethod = 'Z/sqrt(n)';
        if isfield(stats,'zval')
            stat(statcnt).effectSize = stats.zval / ...
                sqrt(sum(cellfun(@numel,ys))); 
        else
            stat(statcnt).effectSize = nan;
        end
        
    end
    
    % n way anova on this statistic
    statcnt = statcnt + 1;
    [pv,tbl,stats] = anovan(yall,{musc,task},'model','interaction',...
                            'varnames',{'muscimol/saline','noise/nonoise'},...
                            'display','off');
    [comp,means,~,gnames] = multcompare(stats,'dimension',[1 2],'display','off');
    stat(statcnt).type = sprintf('%s muscimolxtask anova',stats2run{i});
    stat(statcnt).test = 'anovan';
    stat(statcnt).p = pv;
    stat(statcnt).stats = stats;
    stat(statcnt).table = tbl;
    stat(statcnt).mean = means;
    stat(statcnt).effectSizeMethod = 'partial eta^2';
    for ii = 1:3
        stat(statcnt).effectSize(ii) = tbl{ii+1,2} / (tbl{ii+1,2} + tbl{end-1,2});
    end
    stat(statcnt).multcomp = comp;
    stat(statcnt).multcompGrps = gnames;
    
    set(gca,'xtick',[1 2 4 5]);
    set(gca,'xticklabel',{'Muscimol','Saline','Muscimol','Saline'});
    ylabel(stats2run{i},'Interpreter','none'); plotPrefs;
    xtickangle(45);
    if j > 2
        ylim([0 1]);
    end
    
end

% 2x2x7 ANOVA: intervention (musc/saline) x task (noise/silence) x volume
statcnt = statcnt + 1;
data = T.rate;
intervention = repmat(contains(T.condition,'muscimol'),1,7);
task = repmat(T.control,1,7);
volume = repmat([0 1 2 3 4 5 6],size(data,1),1);
[pv,tbl,stats] = anovan(data(:),{intervention(:),task(:),volume(:)},...
                        'varnames',{'intervention','task','volume'},...
                        'display','on','model','interaction');
[comp,means,~,gnames] = multcompare(stats,'dimension',[1 2],'display','off');
stat(statcnt).type = '3 way ANOVA: musc/saline x noise/silence x volume ~ hit rate';
stat(statcnt).test = 'anovan';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).table = tbl;
stat(statcnt).mean = means;
stat(statcnt).effectSizeMethod = 'partial eta^2';
for i = 1:6
    stat(statcnt).effectSize(i) = tbl{i+1,2} / (tbl{i+1,2} + tbl{end-1,2});
end
for i = 1:2
    [comp,means,~,gnames] = multcompare(stats,'dimension',i,'display','off');
    stat(statcnt).multcomp(i,:) = comp;
    stat(statcnt).multcompGrps{i} = gnames;
end

saveFigPDF(f4,[400 600],'./_plots/_muscimol_control_summary.pdf'); 


if false
    groupVars = {'mouse','condition'};
    T.cumLickP = cumsum(T.plick,2) ./ sum(T.plick,2);
    m = grpstats(T,groupVars,'mean','datavars',statVars);

    f2 = figure(2); clf; hold on; clear pp;
    plot(-90,m.mean_rate(contains(m.condition,'saline'),1),'ok','MarkerFaceColor','k');
    plot(-90,m.mean_rate(contains(m.condition,'muscimol'),1),'ok');
    pp(1) = plot(m.mean_snr(1,:),m.mean_rate(contains(m.condition,'saline'),2:end),...
                 'o-k','MarkerFaceColor','k');
    pp(2) = plot(m.mean_snr(1,:),m.mean_rate(contains(m.condition,'muscimol'),2:end),...
                 'o--k');
    legend(pp,'saline','muscimol')
    xlabel('Target Attenuation (dB rel. 25 dB SNR)');
    ylabel('Response Rate'); plotPrefs;

    saveFigPDF(f2,[300 300],'./_plots/_muscimol_control.pdf');
end





































if false

    %% stats across all mice
    f1 = figure(1); clf; clear p;
    saveFigPDF(f1,[400 950],'./_plots/_figures/muscimol_mice_avg.pdf');


    % for each contrast
    for j = 1:2
        mfill = {'none',cclite(j,:)};
        subplot(5,2,j); hold on;
        for i = 1:length(conds)
            
            I = contains(m.condition,conds{i})' & m.contrast == (j-1);
            x = m.mean_snr(I,:); xf = linspace(min(x(:)),max(x(:)),100);
            y = m.mean_rate(I,2:end);
            [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
            
            plot(mean(x,1)',y','color',cclite(j,:),'linestyle',ls{i}),
            errorbar(mean(x,1),mean(y,1),sem(y,1),'o','color',cclite(j,:),...
                     'markerfacecolor',mfill{i},'markersize',ms);
            p(i) =  plot(xf,mdl(prms,xf),'color',colors{j},'linestyle',ls{i},...
                         'linewidth',1);
            errorbar(min(x(:)) - 5,mean(m.mean_rate(I,1),1),sem(m.mean_rate(I,1),1),...
                     'o','color',cclite(j,:),'markerfacecolor',...
                     mfill{i},'markersize',ms);
            
        end
        legend(p,conds,'location','nw');
        set(gca,'xtick',[-5 mean(x)])
        set(gca,'xticklabels',num2str([-inf mean(x)]'))
        xlabel('Target Volume (dB SNR)'); xlim([-10 30]);
        ylabel('P(respond)'); ylim([0 1]);
        plotPrefs;

    end

    stats2run = {'thresh_rate','slope_rate','mean_maxrate','mean_fa'};
    clear s;
    for k = 1:length(stats2run)
        
        hold on;
        cnt = 1;
        for i = 1:2
            
            subplot(5,2,3+(k-1)); hold on;
            xs = []; ys = [];
            
            for j = 1:2
                
                I = contains(m.condition,conds{j})' & m.contrast == (i-1);            
                x = cnt; y = mean(m.(stats2run{k})(I));
                s(j) = scatter(repmat(cnt,sum(I),1),m.(stats2run{k})(I),ms,'k');
                s(j).MarkerFaceColor = 'w'; s(j).MarkerEdgeColor = 'k';
                xs = [xs repmat(cnt,sum(I),1)];
                ys = [ys m.(stats2run{k})(I)];

                b(cnt) = bar(x,y);
                b(cnt).FaceColor = cc(i,:); b(cnt).FaceAlpha = alphas(j);
                b(cnt).LineStyle = ls{j}; b(cnt).EdgeColor = cc(i,:);
                b(cnt).LineWidth = 1;
                
                cnt = cnt + 1;
                if i == 1 & j == 2
                    cnt = cnt + 1;
                end
                

            end
            
            % plot lines and pval
            plot(xs',ys','k');
            [~,pv] = ttest(ys(:,1),ys(:,2));
            plotPvalSmart(pv,mean(xs),max(ys(:)));
            
        end
        set(gca,'xtick',[1 2 4 5]);
        set(gca,'xticklabel',{'Muscimol','Saline','Muscimol','Saline'});
        ylabel(stats2run{k},'Interpreter','none'); plotPrefs;
        xtickangle(45);
        if k > 2
            ylim([0 1]);
        end
    end


    % for each contrast
    for j = 1:2
        mfill = {'none',cclite(j,:)};
        subplot(5,2,6+j); hold on;
        plot([-5 30],[.5 .5],'k');
        for i = 1:length(conds)
            
            I = contains(m.condition,conds{i})' & m.contrast == (j-1);
            x = m.mean_snr(I,:); xf = linspace(min(x(:)),max(x(:)),100);
            y = m.mean_pc(I,:);
            [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
            
            plot(mean(x,1)',y','color',cclite(j,:),'linestyle',ls{i}),
            errorbar(mean(x,1),mean(y,1),sem(y,1),'o','color',cclite(j,:),...
                     'markerfacecolor',mfill{i},'markersize',ms);
            p(i) = plot(xf,mdl(prms,xf),'color',colors{j},'linestyle',ls{i},...
                        'linewidth',1);
            
        end
        set(gca,'xtick',[mean(x)])
        set(gca,'xticklabels',num2str([mean(x)]'))
        xlabel('Target Volume (dB SNR)'); xlim([-5 30]);
        ylabel('PC_{behavior}'); ylim([.4 1]);
        plotPrefs;

    end


    %% actx data
    % conditions for each mouse
    mList = {'AW142','AW146'};
    wd = '~/chris-lab/projects/gain_behavior';
    datDir = [wd '/_data'];
    resDir = [wd '/_data/muscimol'];
    condLong = {'0.25mg/mL muscimol',...
                'saline'};
    condShort = {'musc25','saline'}; 

    % plot params
    ls = {'-','--'};
    ms = 6;
    mstyle = {'o','o'};
    colors = {'b','r'};
    snrs = mean(m.mean_snr);
    xtick = [snrs];
    xtickl = num2str([snrs]');

    ax(1) = subplot(5,2,9);
    ax(2) = subplot(5,2,10);

    % compare effects over different sessions
    order = [2 1];
    for i = 1:length(order)
        
        I = order(i);
        
        % load data
        fn = fullfile(resDir,...
                      sprintf('%s_%s.mat',mList{I},condShort{I}));
        load(fn)
        
        % for each cell
        clear targMean;
        for j = 1:length(res)
            % extract post session
            targMean(j,:,:) = res(j).auc;
        end
        
        % plot each contrast
        clear x y;
        for j = 1:2
            axes(ax(j)); hold on;
            if i == 1
                plot([-5 30],[.5 .5],'k');
            end
            
            % fit
            x = xtick;
            y = mean(squeeze(targMean(:,j,:)));
            [prms,mdl,thresh] = fitLogGrid(x',y');
            
            % plot
            h = errorBars(snrs,squeeze(targMean(:,j,:)),cclite(j,:));
            h.Marker = mstyle{I}; h.LineStyle = 'none'; h.MarkerSize = ms;
            if I == 2
                h.MarkerFaceColor = cclite(j,:);
            end
            plot(xf,mdl(prms,xf),colors{j},'linewidth',1,'linestyle',ls{i})
            axis tight; ylim([.4 1]);
            
        end
        
    end

    for i = 1:length(ax)
        set(ax(i),'xtick',xtick,'xticklabels',xtickl);
        axes(ax(i)); xlabel('Target Volume (dB SNR)'); ylabel('PC_{actx}');
        plotPrefs;
    end

    saveFigPDF(f1,sz1,'./_plots/_figures/muscimol_mice_avg.pdf');
end


