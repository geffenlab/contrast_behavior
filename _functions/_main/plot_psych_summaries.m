function res = plot_psych_summaries(r,ops)

clear res;

%% plot stuff
nrows = 5; ncols = 6;
col = {'b','r'};
col_rgb = [0 0 1; 1 0 0];
col_lite = {[.7 .7 1],[1 .7 .7]};
xl = [-7.5 27.5];
labels = {'Behavior','Population SVM',...
          'Population Criterion Classifier',...
          'Population AUC','Single Neuron Crit. Classifier',...
          'Single Neuron AUC'};

% x values
x_n = mean(r.vols);
x_nn = mean(r.vols_nn);


%% post-processing
fields = {'beh_rate_adj_PC','svm_rate_adj_PC','critp_adj_PC','auc',...
          'mean_critp_adj_PC','mean_auc'};

% all thresholds outliers
thresh_all = [];
for i = 1:length(fields)
    thresh_all = [thresh_all r.([fields{i} '_fit']).threshold];
end
% sessions with significant population auc
sigAUC = sum(r.auc_sig,2,'omitnan');
% exclude CA121 in high contrast (didn't learn)
badMouse = contains(r.mouse,'CA121') & r.contrastI == 1;
% exclude sessions where highest percent correct is below the
% training criterion
badSess = max(r.beh_rate_adj_PC,[],2,'omitnan') < .8;
include = ~badMouse;





%% AVERAGES

% create figures
f2 = figure(121212); clf;
sz = [1200 1100];
set(f2,'Position',[0 0 sz]);



% precompute behavioral psych curves
subplot(nrows,ncols,1); hold on;
plot(xl,[.5 .5],'k--');
for i = 1:2
    I = r.contrastI==(i-1) & include;
    x_b{i} = r.vols_nn(I,:);
    y_b{i} = r.beh_rate_adj_PC(I,:);
    xf = linspace(min(x_b{i}(:)),max(x_b{i}(:)),100);
    [prms_b(i,:),mdl,thresh_b(i)] = fitLogGrid(x_b{i}(:),y_b{i}(:));
    errorBars(x_nn,r.beh_rate_adj_PC(I,:),col{i},[],[],[],[],...
              'o','LineWidth',1,'MarkerSize',3,'MarkerFaceColor',col{i});
    plot(xf,mdl(prms_b(i,:),xf),'color',col{i},'linewidth',1)
    plot([thresh_b(i) thresh_b(i)],...
         [.5 mdl(prms_b(i,:),thresh_b(i))],...
         'color',col{i},'linewidth',.5)
end
xlim(xl); ylim([.45 1]); plotPrefs; 
title(labels{1}); ylabel('Percent Correct');
xlabel('Target Volume (dB SNR)');

% for the remaining fields
for j = 2:length(fields)
    
    subplot(nrows,ncols,j); hold on;
    plot(xl,[.5 .5],'k--');

    % per contrast
    for i = 1:2
        
        I = r.contrastI==(i-1) & include;
        
        % replot behavior
        errorBars(x_nn,r.beh_rate_adj_PC(I,:),col_lite{i},[],[],[],[],...
                  'o','LineWidth',1,'MarkerSize',3, ...
                  'MarkerFaceColor',col_lite{i});
        plot(xf,mdl(prms_b(i,:),xf),'color',col_lite{i},'linewidth',1)
        plot([thresh_b(i) thresh_b(i)],[.5 mdl(prms_b(i,:),thresh_b(i))],...
             'color',col_lite{i},'linewidth',.5)
        
        % fit neural data
        x = r.vols_nn(I,:);
        y = r.(fields{j})(I,:);
        [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
        errorBars(x_nn,r.(fields{j})(I,:),col{i},[],[],[],[],...
                  'o','LineWidth',1,'MarkerSize',3, ...
                  'MarkerFaceColor',col{i});
        plot(xf,mdl(prms,xf),'color',col{i},'linewidth',1)
        plot([thresh thresh],[.5 mdl(prms,thresh)],'color',col{i},'linewidth',.5)
        
    end
    xlim(xl); ylim([.45 1]); plotPrefs; 
    title(labels{j},'interpreter','none'); drawnow;

end

% plot proportion significant neurons
subplot(nrows,ncols,7); hold on;
for i = 1:2
    
    I = r.contrastI==(i-1) & include;
    x = r.vols_nn(I,:);
    y = r.prop_sig_neurons(I,:);
    [prms,mdl,thresh] = fitLogGrid(x(:),y(:));
    errorBars(x_nn,y,col{i},[],[],[],[],...
              'o','LineWidth',1,'MarkerSize',3, ...
              'MarkerFaceColor',col{i});
    plot(xf,mdl(prms,xf),'color',col{i},'linewidth',1)
    plot([thresh thresh],[0 mdl(prms,thresh)],'color',col{i},'linewidth',.5)
end
ylim([0 1]); xlim(xl); xlabel('Target Volume (dB SNR)');
ylabel('% Significant Neurons');

% get sessions with matched volumes
for i = 1:length(r.vols_nn)
    % check if length is good
    ind = ~isnan(r.beh_rate_adj(i,2:end));
    if sum(ind,2) == 6
        matchSess(i,1) = all(r.vols_nn(i,ind) == [0 5 10 15 20 25],2);
        matchRange(i,1) = range(r.vols_nn(i,ind)) == 25;
    else
        matchSess(i,1) = false;
        matchRange(i,1) = false;

    end
    sessRange(i,1) = range(r.vols_nn(i,ind));
end





%% THRESHOLDS
% plot behavior thresholds against neural, for each mouse
grpvar = 'mouse';

% exclude CA121 in high contrast, include sessions with more than 3
% significant population responses to different volumes, include
% sessions with false alarm rates < .3
% include:  
include = ~badMouse & sum(r.auc_sig,2,'omitnan')>3 & r.beh_rate_adj(:,1) < .3;
stat = 'mean';
%include = ~badMouse & ...
%          sum(r.auc_sig,2,'omitnan')>3 & ...
%          r.beh_rate_adj(:,1) < .3;
% stat = 'median'; %'mean';
% ~any(thresh_all < -10,2) & ...
% grpvar = 'sessionID';

statcnt = 0;
for k = 2:length(fields)
    
    contr = r.contrastI(include);
    beh = r.beh_rate_adj_PC_fit.threshold(include);
    neural = r.([fields{k} '_fit']).threshold(include);
    grps = r.(grpvar)(include);
    
    ci = grpstats(contr,{grps,contr},stat)+1;
    cv = col_rgb(ci,:);
    x = grpstats(neural,{grps,contr},stat);
    [y,g] = grpstats(beh,{grps,contr},{stat,'gname'});
    g = g(:,1);
    
    figure(f2);
    subplot(nrows,ncols,6+k); hold on
    plot([0 20],[0 20],'k');

    
    %%
    % linear model for all data
    [b,pv,xp,yp,ypci,lm] = fitlmWrapper(x,y);
    plotlmWrapper(xp,yp,ypci);
    sh = scatter(x,y,20,cv);
    sh.MarkerFaceColor = 'flat';
    sh.MarkerEdgeColor = 'w';
    xlim([0 25]); ylim([0 25]);
    set(gca,'xtick',[0 10 20]);
    set(gca,'ytick',[0 10 20]);
    ylabel('Behavioral Threshold (dB SNR)');
    xlabel('Neural Threshold (dB SNR)');
    plotPrefs; axis tight;
    
    % save stats for linear model
    res.(fields{k}).stats(1).type = 'threshold linear model all data points';
    res.(fields{k}).stats(1).test = 'fitlm';
    res.(fields{k}).stats(1).p = pv;
    res.(fields{k}).stats(1).stats = lm;
    res.(fields{k}).stats(1).n = numel(x);
    

    %%
    % n way anova for unmatched data
    data = [x;y];
    measure = [ones(size(x)); ones(size(y))*2];
    contrast = [ci;ci];
    [pv3 tbl stats] = anovan(data,{measure,contrast},...
                             'model','interaction',...
                             'varnames',{'threshold measure','contrast'},...
                             'display','off');
    res.(fields{k}).stats(2).type = 'two way anova of contrast x threshold';
    res.(fields{k}).stats(2).test = 'anovan';
    res.(fields{k}).stats(2).p = pv3;
    res.(fields{k}).stats(2).stats = tbl;
    res.(fields{k}).stats(2).n = numel(x);
    res.(fields{k}).stats(2).eta2.threshold_measure = tbl{2,2} / tbl{6,2};
    res.(fields{k}).stats(2).eta2.contrast = tbl{3,2} / tbl{6,2};
    res.(fields{k}).stats(2).eta2.interaction = tbl{4,2} / tbl{6,2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1],'display','off');
    res.(fields{k}).stats(2).multcomp(1).comps = m;
    res.(fields{k}).stats(2).multcomp(1).means = c;
    res.(fields{k}).stats(2).multcomp(1).gnames = gnames;
    [c,m,h,gnames] = multcompare(stats,'dimension',[2],'display','off');
    res.(fields{k}).stats(2).multcomp(2).comps = m;
    res.(fields{k}).stats(2).multcomp(2).means = c;
    res.(fields{k}).stats(2).multcomp(2).gnames = gnames;
    
    %%
    % mixed effects model linear model with contrast predictor
    t = table();
    t.neural = x;
    t.behavior = y;
    t.contrast = ci;
    t.mouse = g;
    lme = fitlme(t,'behavior ~ neural + contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_neural = fitlme(t,'behavior ~ contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_contrast = fitlme(t,'behavior ~ neural + (contrast-1|mouse)',...
                 'StartMethod','random');
    
    res.(fields{k}).stats(3).type = ['beh_thresh = neural_thresh + ' ...
                        'contrast + (contrast-1|mouse)'];
    res.(fields{k}).stats(3).test = 'fitlme';
    res.(fields{k}).stats(3).stats = lme;
    res.(fields{k}).stats(3).coefficients = lme.Coefficients;
    res.(fields{k}).stats(3).n = numel(x);
    res.(fields{k}).stats(3).mdlcomp_no_neural = ...
        compare(lme_no_neural,lme);
    res.(fields{k}).stats(3).mdlcomp_no_contrast = ...
        compare(lme_no_contrast,lme);
    
    if (2 + 2) == 5
        X = [x ci];
        lm = fitlm(X,y,'linear');
        res.(fields{k}).stats(3).type = 'beh_thresh = neural_thresh + contrast';
        res.(fields{k}).stats(3).test = 'fitlm';
        res.(fields{k}).stats(3).stats = lm;
        res.(fields{k}).stats(3).n = numel(x);
    end
    
    title(sprintf('p_{neural} = %04.3f\np_{contrast} = %04.3f',...
                  res.(fields{k}).stats(3).mdlcomp_no_neural.pValue(2),...
                  res.(fields{k}).stats(3).mdlcomp_no_contrast.pValue(2)));
    drawnow;
    
end



%% SLOPES
lims = [.025 .15];
for k = 2:length(fields)
    
    contr = r.contrastI(include);
    beh = r.beh_rate_adj_PC_fit.max_slope(include);
    neural = r.([fields{k} '_fit']).max_slope(include);
    grps = r.(grpvar)(include);
    
    ci = grpstats(contr,{grps,contr},'mean')+1;
    cv = col_rgb(ci,:);
    x = grpstats(neural,{grps,contr},'mean');
    [y,g] = grpstats(beh,{grps,contr},{'mean','gname'});
    g = g(:,1);
    
    subplot(nrows,ncols,12+k); hold on
    plot(lims,lims,'k');

    
    %%
    % linear model for all data
    [b,pv,xp,yp,ypci,lm] = fitlmWrapper(x,y);
    plotlmWrapper(xp,yp,ypci);
    sh = scatter(x,y,20,cv);
    sh.MarkerFaceColor = 'flat';
    sh.MarkerEdgeColor = 'w';
    xlim(lims); ylim(lims);
    %set(gca,'xtick',[0 10 20]);
    %set(gca,'ytick',[0 10 20]);
    ylabel('Behavioral Slope (PC/dB)');
    xlabel('Neural Slope (PC/dB)');
    plotPrefs; axis tight;
    
    % save stats for linear model
    res.(fields{k}).stats(4).type = 'slope linear model all data points';
    res.(fields{k}).stats(4).test = 'fitlm';
    res.(fields{k}).stats(4).p = pv;
    res.(fields{k}).stats(4).stats = lm;
    res.(fields{k}).stats(4).n = numel(x);
    
    
    %%
    % n-way anova
    data = [x;y];
    measure = [ones(size(x)); ones(size(y))*2];
    contrast = [ci;ci];
    [pv3 tbl stats] = anovan(data,{measure,contrast},...
                             'model','interaction',...
                             'varnames',{'slope measure','contrast'},...
                             'display','off');
    res.(fields{k}).stats(5).type = 'two way anova of contrast x slope';
    res.(fields{k}).stats(5).test = 'anovan';
    res.(fields{k}).stats(5).p = pv3;
    res.(fields{k}).stats(5).stats = tbl; 
    res.(fields{k}).stats(5).n = numel(x);
    res.(fields{k}).stats(5).eta2.threshold_measure = tbl{2,2} / tbl{6,2};
    res.(fields{k}).stats(5).eta2.contrast = tbl{3,2} / tbl{6,2};
    res.(fields{k}).stats(5).eta2.interaction = tbl{4,2} / tbl{6,2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1],'display','off');
    res.(fields{k}).stats(5).multcomp(1).comps = m;
    res.(fields{k}).stats(5).multcomp(1).means = c;
    res.(fields{k}).stats(5).multcomp(1).gnames = gnames;
    [c,m,h,gnames] = multcompare(stats,'dimension',[2],'display','off');
    res.(fields{k}).stats(5).multcomp(2).comps = m;
    res.(fields{k}).stats(5).multcomp(2).means = c;
    res.(fields{k}).stats(5).multcomp(2).gnames = gnames;
    
    
    %%
    % mixed effects model linear model with contrast predictor
    t = table();
    t.neural = x;
    t.behavior = y;
    t.contrast = ci;
    t.mouse = g;
    lme = fitlme(t,'behavior ~ neural + contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_neural = fitlme(t,'behavior ~ contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_contrast = fitlme(t,'behavior ~ neural + (contrast-1|mouse)',...
                 'StartMethod','random');
    
    res.(fields{k}).stats(6).type = ['beh_slope = neural_slope + ' ...
                        'contrast + (contrast-1|mouse)'];
    res.(fields{k}).stats(6).test = 'fitlme';
    res.(fields{k}).stats(6).stats = lme;
    res.(fields{k}).stats(6).coefficients = lme.Coefficients;
    res.(fields{k}).stats(6).n = numel(x);
    res.(fields{k}).stats(6).mdlcomp_no_neural = ...
        compare(lme_no_neural,lme);
    res.(fields{k}).stats(6).mdlcomp_no_contrast = ...
        compare(lme_no_contrast,lme);
    
    if (2 + 2) == 5
        X = [x ci];
        lm = fitlm(X,y,'linear');
        res.(fields{k}).stats(6).type = 'beh_slope = neural_slope + contrast';
        res.(fields{k}).stats(6).test = 'fitlm';
        res.(fields{k}).stats(6).stats = lm;
        res.(fields{k}).stats(6).n = numel(x);
    end
    
    title(sprintf('p_{neural} = %04.3f\np_{contrast} = %04.3f',...
                  res.(fields{k}).stats(6).mdlcomp_no_neural.pValue(2),...
                  res.(fields{k}).stats(6).mdlcomp_no_contrast.pValue(2)));
    drawnow;
    
    
    %% n-way anova only for mice with matched targets
    % set up mice
    ind = matchRange;
    contr = r.contrastI(ind);
    beh = r.beh_rate_adj_PC_fit.max_slope(ind);
    neural = r.([fields{k} '_fit']).max_slope(ind);
    grps = r.(grpvar)(ind);
    
    % anova
    ci = grpstats(contr,{grps,contr},'mean')+1;
    cv = col_rgb(ci,:);
    x = grpstats(neural,{grps,contr},'mean');
    [y,g,nn] = grpstats(beh,{grps,contr},{'mean','gname','numel'});
    data = [x;y];
    measure = [ones(size(x)); ones(size(y))*2];
    contrast = [ci;ci];
    [pv3 tbl stats] = anovan(data,{measure,contrast},...
                             'model','interaction',...
                             'varnames',{'slope measure','contrast'},...
                             'display','off');    
    res.(fields{k}).stats(7).type = ['two way anova of contrast x ' ...
                        'slope with matched volumes'];
    res.(fields{k}).stats(7).test = 'anovan';
    res.(fields{k}).stats(7).p = pv3;
    res.(fields{k}).stats(7).stats = tbl; 
    res.(fields{k}).stats(7).n = numel(x);
    res.(fields{k}).stats(7).eta2.threshold_measure = tbl{2,2} / tbl{6,2};
    res.(fields{k}).stats(7).eta2.contrast = tbl{3,2} / tbl{6,2};
    res.(fields{k}).stats(7).eta2.interaction = tbl{4,2} / tbl{6,2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1],'display','off');
    res.(fields{k}).stats(7).multcomp(1).comps = m;
    res.(fields{k}).stats(7).multcomp(1).means = c;
    res.(fields{k}).stats(7).multcomp(1).gnames = gnames;
    [c,m,h,gnames] = multcompare(stats,'dimension',[2],'display','off');
    res.(fields{k}).stats(7).multcomp(2).comps = m;
    res.(fields{k}).stats(7).multcomp(2).means = c;
    res.(fields{k}).stats(7).multcomp(2).gnames = gnames;
    
    
    % lmes for matched mice
    t = table();
    t.neural = x;
    t.behavior = y;
    t.contrast = ci;
    t.mouse = cat(1,g(:,1));
    lme = fitlme(t,'behavior ~ neural + contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_neural = fitlme(t,'behavior ~ contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_contrast = fitlme(t,'behavior ~ neural + (contrast-1|mouse)',...
                 'StartMethod','random');
    
    res.(fields{k}).stats(11).type = ['beh_slope = neural_slope + ' ...
                        'contrast + (contrast-1|mouse) *MATCHED TARGS*'];
    res.(fields{k}).stats(11).test = 'fitlme';
    res.(fields{k}).stats(11).stats = lme;
    res.(fields{k}).stats(11).coefficients = lme.Coefficients;
    res.(fields{k}).stats(11).n = numel(x);
    res.(fields{k}).stats(11).mdlcomp_no_neural = ...
        compare(lme_no_neural,lme);
    res.(fields{k}).stats(11).mdlcomp_no_contrast = ...
        compare(lme_no_contrast,lme);
    
    if (2 + 2) == 5
        X = [x ci];
        lm = fitlm(X,y,'linear');
        res.(fields{k}).stats(11).type = 'beh_slope = neural_slope + contrast';
        res.(fields{k}).stats(11).test = 'fitlm';
        res.(fields{k}).stats(11).stats = lm;
        res.(fields{k}).stats(11).n = numel(x);
    end
    
    drawnow;
    
        
    % plot
    subplot(nrows,ncols,24+k); hold on
    plot([.01 .055],[.01 .055],'k');
    [b,pv,xp,yp,ypci,lm] = fitlmWrapper(x,y);
    plotlmWrapper(xp,yp,ypci);
    sh = scatter(x,y,20,cv);
    sh.MarkerFaceColor = 'flat';
    sh.MarkerEdgeColor = 'w';
    xlim(lims); ylim(lims);
    %set(gca,'xtick',[0 10 20]);
    %set(gca,'ytick',[0 10 20]);
    ylabel('Behavioral Slope (PC/dB)');
    xlabel('Neural Slope (PC/dB)');
    plotPrefs; axis tight;
    
    title(sprintf('p_{neural} = %04.3f\np_{contrast} = %04.3f',...
                  res.(fields{k}).stats(11).mdlcomp_no_neural.pValue(2),...
                  res.(fields{k}).stats(11).mdlcomp_no_contrast.pValue(2)));
    
end




%% SENSITIVITY
% plot behavior sensitivity against neural, for each mouse
lims = [0 3];
for k = 2:length(fields)
    
    contr = r.contrastI(include);
    beh = r.beh_rate_adj_PC_fit.sensitivity(include);
    neural = r.([fields{k} '_fit']).sensitivity(include);
    grps = r.(grpvar)(include);
    
    ci = grpstats(contr,{grps,contr},'mean')+1;
    cv = col_rgb(ci,:);
    x = grpstats(neural,{grps,contr},'mean');
    [y,g] = grpstats(beh,{grps,contr},{'mean','gname'});
    g = g(:,1);
    
    subplot(nrows,ncols,18+k); hold on
    plot(lims,lims,'k');

    
    %%
    % linear model for all data
    [b,pv,xp,yp,ypci,lm] = fitlmWrapper(x,y);
    plotlmWrapper(xp,yp,ypci);
    sh = scatter(x,y,20,cv);
    sh.MarkerFaceColor = 'flat';
    sh.MarkerEdgeColor = 'w';
    xlim(lims); ylim(lims);
    %set(gca,'xtick',[0 10 20]);
    %set(gca,'ytick',[0 10 20]);
    ylabel('Behavioral Sensitivity');
    xlabel('Neural Sensitivity');
    plotPrefs; axis tight;
    
    % save stats for linear model
    res.(fields{k}).stats(8).type = 'sensitivity linear model all data points';
    res.(fields{k}).stats(8).test = 'fitlm';
    res.(fields{k}).stats(8).p = pv;
    res.(fields{k}).stats(8).stats = lm;
    res.(fields{k}).stats(8).n = numel(x);
    
    
    %%
    % format data for rmanova
    data = [x;y];
    measure = [ones(size(x)); ones(size(y))*2];
    contrast = [ci;ci];
    [pv3 tbl stats] = anovan(data,{measure,contrast},...
                             'model','interaction',...
                             'varnames',{'threshold measure','contrast'},...
                             'display','off');
    res.(fields{k}).stats(9).type = 'two way anova of contrast x sensitivity';
    res.(fields{k}).stats(9).test = 'anovan';
    res.(fields{k}).stats(9).p = pv3;
    res.(fields{k}).stats(9).stats = stats;
    res.(fields{k}).stats(9).table = tbl;
    res.(fields{k}).stats(9).n = numel(x);
    res.(fields{k}).stats(9).eta2.threshold_measure = tbl{2,2} / tbl{6,2};
    res.(fields{k}).stats(9).eta2.contrast = tbl{3,2} / tbl{6,2};
    res.(fields{k}).stats(9).eta2.interaction = tbl{4,2} / tbl{6, ...
                        2};
    [c,m,h,gnames] = multcompare(stats,'dimension',[1],'display','off');
    res.(fields{k}).stats(9).multcomp(1).comps = m;
    res.(fields{k}).stats(9).multcomp(1).means = c;
    res.(fields{k}).stats(9).multcomp(1).gnames = gnames;
    [c,m,h,gnames] = multcompare(stats,'dimension',[2],'display','off');
    res.(fields{k}).stats(9).multcomp(2).comps = m;
    res.(fields{k}).stats(9).multcomp(2).means = c;
    res.(fields{k}).stats(9).multcomp(2).gnames = gnames;
    
    %%
    % mixed effects model linear model with contrast predictor
    t = table();
    t.neural = x;
    t.behavior = y;
    t.contrast = ci;
    t.mouse = g;
    lme = fitlme(t,'behavior ~ neural + contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_neural = fitlme(t,'behavior ~ contrast + (contrast-1|mouse)',...
                 'StartMethod','random');
    lme_no_contrast = fitlme(t,'behavior ~ neural + (contrast-1|mouse)',...
                 'StartMethod','random');
    
    res.(fields{k}).stats(10).type = ['beh_sense = neural_sense + ' ...
                        'contrast + (contrast-1|mouse)'];
    res.(fields{k}).stats(10).test = 'fitlme';
    res.(fields{k}).stats(10).stats = lme;
    res.(fields{k}).stats(10).coefficients = lme.Coefficients;
    res.(fields{k}).stats(10).n = numel(x);
    res.(fields{k}).stats(10).mdlcomp_no_neural = ...
        compare(lme_no_neural,lme);
    res.(fields{k}).stats(10).mdlcomp_no_contrast = ...
        compare(lme_no_contrast,lme);
    
    if (2+2) == 5
        % linear model with contrast predictor
        X = [x ci];
        lm = fitlm(X,y,'linear');
        res.(fields{k}).stats(10).type = 'beh_sense = neural_sense + contrast';
        res.(fields{k}).stats(10).test = 'fitlm';
        res.(fields{k}).stats(10).stats = lm;
        res.(fields{k}).stats(10).n = numel(x);
    end
    
    title(sprintf('p_{neural} = %04.3f\np_{contrast} = %04.3f',...
                  res.(fields{k}).stats(10).mdlcomp_no_neural.pValue(2),...
                  res.(fields{k}).stats(10).mdlcomp_no_contrast.pValue(2)));
    
end

saveFigPDF(f2,sz,'./_plots/_psych_summary.pdf')






%%
% average curve for each mouse for each stat
%clear res;
uM = unique(r.mouse);
for k = 1:length(fields)
    for i = 1:length(uM)
        for j = 1:2
            
            fprintf('%s: %d %d\n',fields{k},i,j);
            
            I = r.contrastI==(j-1) & contains(r.mouse,uM{i}) & ...
                include;
            
            
            res.(fields{k}).params(i,j,:) = nan(1,4);
            res.(fields{k}).thresh(i,j) = nan;
            res.(fields{k}).thresh75(i,j) = nan;
            res.(fields{k}).sensitivity(i,j) = nan;
            res.(fields{k}).maxslope(i,j) = nan;

            if sum(I) > 0
                
                % fit curve
                x = r.vols_nn(I,:);
                y = r.(fields{k})(I,:);
                xn = mean(x,1,'omitnan');  yn = mean(y,1,'omitnan');
                xn = xn(~isnan(yn)); yn = yn(~isnan(yn));
                [prms,mdl,thresh,sense,~,~,thresh75] = fitLogGrid(xn,yn,[],[],[],.75);
                
                % max slope
                mxslope = max(diff(yn)./diff(xn));
                
                % results
                res.(fields{k}).params(i,j,:) = prms;
                res.(fields{k}).thresh(i,j) = thresh;
                res.(fields{k}).thresh75(i,j) = thresh75;
                res.(fields{k}).sensitivity(i,j) = sense;
                res.(fields{k}).maxslope(i,j) = mxslope;
                res.(fields{k}).x{i,j} = x;
                res.(fields{k}).y{i,j} = y;

                
            end
            
        end
    end
end



%% behavior plots of individual mice
f1 = figure(12412); clf;

for i = 1:length(uM)
    for k = 1:length(fields)
        
        subplot(length(uM),length(fields),...
                (i-1)*(length(fields))+k);
        plot(xl,[.5 .5],'k--');
        hold on;
        
        for j = 1:2
            
            
            I = r.contrastI==(j-1) & contains(r.mouse,uM{i}) & ...
                include;
            
            if sum(I) > 0
                if k == 1
                    % plot only behavior
                    x = r.vols_nn(I,:);
                    y = r.(fields{1})(I,:);
                    xn = mean(x,1,'omitnan');  yn = mean(y,1,'omitnan');
                    xn = xn(~isnan(yn)); yn = yn(~isnan(yn));
                    xf = linspace(min(xn(:)),max(xn(:)),100);
                    
                    scatter(xn,yn,50,col_lite{j},'.')
                    plot(xf,mdl(res.(fields{1}).params(i,j,:),xf), ...
                         col{j},'linewidth',1);
                    thr = res.(fields{1}).thresh(i,j);
                    plot([thr,thr],[.5 mdl(res.(fields{1}).params(i,j,:),thr)],...
                         col{j},'linewidth',.5);
                    
                else
                    % get performance for this measure
                    x = r.vols_nn(I,:);
                    y = r.(fields{k})(I,:);
                    xn = mean(x,1,'omitnan');  yn = mean(y,1,'omitnan');
                    xn = xn(~isnan(yn)); yn = yn(~isnan(yn));
                    xf = linspace(min(xn(:)),max(xn(:)),100);
                    
                    % plot behavior
                    plot(xf,mdl(res.(fields{1}).params(i,j,:),xf),...
                         'color',col_lite{j},'linewidth',1);
                    thr = res.(fields{1}).thresh(i,j);
                    plot([thr,thr],[.5 mdl(res.(fields{1}).params(i,j,:),thr)],...
                         'color',col_lite{j},'linewidth',.5);
                    
                    % plot measure
                    scatter(xn,yn,50,col_lite{j},'.');
                    plot(xf,mdl(res.(fields{k}).params(i,j,:),xf),...
                         'color',col{j},'linewidth',1);
                    thr = res.(fields{k}).thresh(i,j);
                    plot([thr,thr],[.5 mdl(res.(fields{k}).params(i,j,:),thr)],...
                         'color',col{j},'linewidth',.5);
                    xlim(xl); ylim([.4 1]); plotPrefs;
                    
                    
                end
            end
        end
        
        if k == 1 & i == length(uM)
            xlabel('Target Volume');
        end
        if k == 1
            ylabel(sprintf('%s\nPC',uM{i}));
        end
        if i == 1
            title(sprintf('%s',labels{k}));
        end
        
    end
end


saveFigPDF(f1,[800 1300],'./_plots/_psych_individual_curves.pdf')









