function stat = plot_gain_psych_stats(rln,r_psych,g,include,stat,statcnt,facut)

% group gain by sessions
grp = cat(1,rln.sessionID);
[sess_gain,gain_sess] = grpstats(g(include,:),grp(include),{'mean','gname'});

% find groups in common and align the data to each other
gain_sess = cellfun(@str2num,gain_sess);
psych_sess = r_psych.sessionID;
[C,ia,ib] = intersect(gain_sess,psych_sess);
sess_gain = sess_gain(ia,:);
thresh = r_psych.beh_rate_adj_PC_fit.threshold(ib);
slope = r_psych.beh_rate_adj_PC_fit.max_slope(ib);
psy_mouse = r_psych.mouse(ib);
sessions = gain_sess(ia);

% contrast colors
contrast = r_psych.contrastI(ib);
cols = [0 0 1; 1 0 0];
cv = cols(contrast+1,:);
cva = cols(~contrast+1,:); % flipped for adaptation
nrows = 7;
ncols = 2;

pvols = r_psych.vols_nn(ib,:);
prate = r_psych.beh_rate_adj(ib,2:end);

% get sessions with matched volumes
for i = 1:length(pvols)
    % check if length is good
    ind = ~isnan(prate(i,:));
    if sum(ind,2) == 6
        matchSess(i,1) = all(pvols(i,ind) == [0 5 10 15 20 25],2);
    else
        matchSess(i,1) = false;
    end
end

% what is the mouse for each session?
clear gain_mouse;
for i = 1:length(sessions)
    I = find(cat(1,rln.sessionID) == sessions(i),1,'first');
    str = strsplit(rln(I).cellID,'_');
    gain_mouse{i} = str{1};
end

% table for repeated measures
t = table;
t.mouse = gain_mouse';
t.contrast = contrast;
t.gain_adapt = sess_gain(:,1);
t.gain_target = sess_gain(:,2);
t.thresh = thresh;
t.slope = slope;




%% per session results

% session inclusion
[~,slope_out] = rmoutliers(t.slope);
[~,thresh_out] = rmoutliers(t.thresh);
fa = r_psych.beh_rate_adj(ib,1);
incl = fa < facut;


%
% mixed effects model of threshold using contrast,gain as fixed
% effects, mouse ID as random effect (TARGET)
formula = 'thresh ~ contrast + gain_target + (contrast-1|mouse)';
lme = fitlme(t(incl,:),formula);
lme_contrast = fitlme(t(incl,:),'thresh ~ contrast + (contrast-1|mouse)');
lme_gain = fitlme(t(incl,:),'thresh ~ gain_target + (contrast-1|mouse)');

subplot(nrows,ncols,1);
scatter(t.gain_target(incl),t.thresh(incl),40,cv(incl,:),'.');
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Threshold (dB SNR)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of threshold ~ contrast + target_gain over sessions';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(t,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;
stat(statcnt).lme_contrast = lme_contrast;
stat(statcnt).lme_gain = lme_gain;
stat(statcnt).mdlcomp_contrast = compare(lme_contrast,lme);
stat(statcnt).mdlcomp_gain = compare(lme_gain,lme);


%
% mixed effects model of slope using contrast,gain as fixed
% effects, mouse ID as random effect (TARGET)
formula = 'slope ~ 1 + contrast + gain_target + (contrast-1|mouse)';
lme = fitlme(t(incl,:),formula);
lme_contrast = fitlme(t(incl,:),'slope ~ 1 + contrast + (contrast-1|mouse)');
lme_gain = fitlme(t(incl,:),'slope ~ 1 + gain_target + (contrast-1|mouse)');

subplot(nrows,ncols,2);
scatter(t.gain_target(incl),t.slope(incl),40,cv(incl,:),'.');
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Slope (PC/dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of slope ~ contrast + target_gain over sessions';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(t,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;
stat(statcnt).lme_contrast = lme_contrast;
stat(statcnt).lme_gain = lme_gain;
stat(statcnt).mdlcomp_contrast = compare(lme_contrast,lme);
stat(statcnt).mdlcomp_gain = compare(lme_gain,lme);

% ****uncomment me if you want to look at only matched sessions
%keyboard
%formula = 'slope ~ 1 + contrast + gain_target + (1|mouse)';
%lme = fitlme(t(incl&matchSess,:),formula);
%lme_contrast = fitlme(t(incl&matchSess,:),'slope ~ 1 + contrast + (1|mouse)');
%lme_gain = fitlme(t(incl&matchSess,:),'slope ~ 1 + gain_target + (1|mouse)');
%scatter(t.gain_target(incl&matchSess),t.slope(incl&matchSess),40,cv(incl&matchSess,:),'.');


%
% mixed effects model of threshold using contrast,gain as fixed
% effects, mouse ID as random effect (ADAPTATION)
formula = 'thresh ~ 1 + contrast + gain_adapt + (contrast-1|mouse)';
lme = fitlme(t(incl,:),formula);
lme_contrast = fitlme(t(incl,:),'thresh ~ 1 + contrast + (contrast-1|mouse)');
lme_gain = fitlme(t(incl,:),'thresh ~ 1 + gain_adapt + (contrast-1|mouse)');

subplot(nrows,ncols,5);
scatter(t.gain_adapt(incl),t.thresh(incl),40,cva(incl,:),'.');
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Threshold (dB SNR)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of threshold ~ contrast + adapt_gain over sessions';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(t,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;
stat(statcnt).lme_contrast = lme_contrast;
stat(statcnt).lme_gain = lme_gain;
stat(statcnt).mdlcomp_contrast = compare(lme_contrast,lme);
stat(statcnt).mdlcomp_gain = compare(lme_gain,lme);




%
% mixed effects model of slope using contrast,gain as fixed
% effects, mouse ID as random effect (ADAPTATION)
formula = 'slope ~ 1 + contrast + gain_adapt + (contrast-1|mouse)';
lme = fitlme(t(incl,:),formula);
lme_contrast = fitlme(t(incl,:),'slope ~ 1 + contrast + (contrast-1|mouse)');
lme_gain = fitlme(t(incl,:),'slope ~ 1 + gain_adapt + (contrast-1|mouse)');

subplot(nrows,ncols,6);
scatter(t.gain_adapt(incl),t.slope(incl),40,cva(incl,:),'.');
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Slope (PC/dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of slope ~ contrast + adapt_gain over sessions';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(t,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;
stat(statcnt).lme_contrast = lme_contrast;
stat(statcnt).lme_gain = lme_gain;
stat(statcnt).mdlcomp_contrast = compare(lme_contrast,lme);
stat(statcnt).mdlcomp_gain = compare(lme_gain,lme);



%%
% median split of data based on gain, plot mean psychometric
% functions
vols = r_psych.vols_nn(ib,:);
pc = r_psych.beh_rate_adj_PC(ib,:);

% adapt
prct = prctile(t.gain_adapt(t.contrast == 0 & incl),[0 50 100]);
%clrs = repmat(linspace(.6,0,length(prct)-1)',1,3);
clrs = [.6 .6 1; 0 0 1];
bins = linspace(prct(1),prct(end),20);

subplot(nrows,ncols,9); hold on;
plot([-7.5 27.5],[.5 .5],'k');
for i = 1:length(prct)-1
    I = t.gain_adapt >= prct(i) & ...
        t.gain_adapt < prct(i+1) & ...
        t.contrast == 0 & incl; 
    
    subplot(nrows,ncols,9); hold on;
    x = mean(vols(I,2:end),1,'omitnan'); xf = linspace(min(x),max(x),100);
    y = mean(pc(I,2:end),1,'omitnan');
    eh(i) = errorBars(x,pc(I,2:end),clrs(i,:),'sem',[],[],[],...
              'o','LineWidth',1,'MarkerSize',5, ...
              'MarkerFaceColor',clrs(i,:));
    [prms,mdl,thr] = fitLogGrid(x,y);
    %plot(x,y,'.','color',clrs(i,:));
    plot(xf,mdl(prms,xf),'color',clrs(i,:),'linewidth',1);
    plot([thr thr],[.5 mdl(prms,thr)],...
         '--','color',clrs(i,:),'linewidth',1);
    
    subplot(nrows,ncols,11); hold on;
    histogram(t.gain_adapt(I),bins,'FaceColor',clrs(i,:));
    plotPrefs;
    
    
end
subplot(nrows,ncols,11); hold on;
plot([prct(2) prct(2)],ylim,'r--','linewidth',1);
xlabel('Gain_{adapt}');

subplot(nrows,ncols,9); hold on;
axis tight; ylim([.4 1]);
legend(eh,{'Low Gain','High Gain'},'location','se');
title('Adapt Gain');
plotPrefs;

% stats
I1 = t.gain_adapt >= prct(2) & ...
     t.gain_adapt < prct(3) & ...
        t.contrast == 0 & incl;
g_high = repmat(I1,1,size(pc,2)-1);
y = pc(:,2:end);
x = vols(:,2:end);
m = t.mouse;
m = repmat(m,1,size(pc,2)-1);
[pv,tbl,stats] = anovan(y(:),{x(:),g_high(:),m(:)},'random',3,...
                        'varnames',{'volume','gain','mouse'},...
                        'display','off');
multcompare(stats,'dimension',2,'display','off');



% target
prct = prctile(t.gain_target(t.contrast == 0 & incl),[0 50 100]);
bins = linspace(prct(1),prct(end),20);

subplot(nrows,ncols,10); hold on;
plot([-7.5 27.5],[.5 .5],'k');
for i = 1:length(prct)-1
    
    I = t.gain_target >= prct(i) & ...
        t.gain_target < prct(i+1) & ...
        t.contrast == 0 & incl; 
    
    subplot(nrows,ncols,10); hold on;
    x = mean(vols(I,2:end),1,'omitnan'); xf = linspace(min(x),max(x),100);
    y = mean(pc(I,2:end),1,'omitnan');
    eh(i) = errorBars(x,pc(I,2:end),clrs(i,:),'sem',[],[],[],...
              'o','LineWidth',1,'MarkerSize',5, ...
              'MarkerFaceColor',clrs(i,:));
    [prms,mdl,thr] = fitLogGrid(x,y);
    %plot(x,y,'.','color',clrs(i,:));
    plot(xf,mdl(prms,xf),'color',clrs(i,:),'linewidth',1);
    plot([thr thr],[.5 mdl(prms,thr)],...
         '--','color',clrs(i,:),'linewidth',1);
    
    subplot(nrows,ncols,12); hold on;
    histogram(t.gain_target(I),bins,'FaceColor',clrs(i,:));
    plotPrefs;
    
end
subplot(nrows,ncols,12); hold on;
plot([prct(2) prct(2)],ylim,'r--','linewidth',1);
xlabel('Gain_{target}');

subplot(nrows,ncols,10);
legend(eh,{'Low Gain','High Gain'},'location','se');
axis tight; ylim([.4 1]);
title('Target Gain');
plotPrefs;

% stats
I1 = t.gain_target >= prct(2) & ...
     t.gain_target < prct(3) & ...
        t.contrast == 0 & incl;
g_high = repmat(I1,1,size(pc,2)-1);
y = pc(:,2:end);
x = vols(:,2:end);
m = t.mouse;
m = repmat(m,1,size(pc,2)-1);
[pv,tbl,stats] = anovan(y(:),{x(:),g_high(:),m(:)},'random',3,...
                        'varnames',{'volume','gain','mouse'},...
                        'display','off');
multcompare(stats,'dimension',2,'display','off');



%% target vs adaptation period gain across contrasts

%  gg = [t.gain_adapt; t.gain_target];
%  cc = [t.contrast; t.contrast];
%  prd = [zeros(size(t,1),1); ones(size(t,1),1)];
%  mice = [t.mouse; t.mouse];
%  index = [incl; incl];
%  
%  [pv,tbl,stats] = anovan(gg(index),{cc(index),prd(index)},'model',3,...
%                          'varnames',{'Contrast','Adapt/Target'},...
%                          'display','off');





%% per mouse results
% compute average gain of included sessions
ct_measure = 'median';
if contains(ct_measure,'median');
    ct_measure = @(x) median(x,1);
end
% ct_measure = 'mean';
[mouse_gain,gain_grp,grp_cnt] = grpstats(sess_gain(incl,:),{gain_mouse(incl)',contrast(incl)},...
                                 {ct_measure,'gname','numel'});

% compute threshold for each mouse
[mouse_thresh,thresh_grp] = grpstats(thresh(incl),{psy_mouse(incl)',contrast(incl)},...
                                     {ct_measure,'gname'});

% compute slope for each mouse
[mouse_slope,slope_grp] = grpstats(slope(incl),{psy_mouse(incl)',contrast(incl)},...
                                     {ct_measure,'gname'});

% mouse colors
cc = cellfun(@str2num,gain_grp(:,2));
mcv = cols(cc+1,:);
mcva = cols(~cc+1,:);

% mouse table
m = table;
m.mouse = gain_grp(:,1);
m.contrast = cc;
m.thresh = mouse_thresh;
m.slope = mouse_slope;
m.gain_adapt = mouse_gain(:,1);
m.gain_target = mouse_gain(:,2);

% remove CA121 who didn't perform in high contrast
out = contains(slope_grp(:,2),'1') & contains(slope_grp(:,1),'CA121');
%out = all(grp_cnt==1,2);


%
% regress threshold against contrast and gain (TARGET)
y = m.thresh;
[~,out] = rmoutliers(y);
%formula = 'thresh ~ contrast + gain_target + (1|mouse)';
formula = 'thresh ~ contrast + gain_target';
lme = fitlme(m(~out,:),formula);

subplot(nrows,ncols,3)
xx = m.gain_target(~out); yy = y(~out);
sh = scatter(xx,yy,20,mcv(~out,:));
sh.MarkerFaceColor = 'flat';
sh.MarkerEdgeColor = 'w';
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Threshold (dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of threshold ~ contrast + gain over mice';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(y,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;


%
% regress slope against contrast and gain (TARGET)
y = m.slope;
[~,out] = rmoutliers(y);
%formula = 'slope ~ contrast + gain_target + (1|mouse)';
formula = 'thresh ~ contrast + gain_target';
lme = fitlme(m(~out,:),formula);

subplot(nrows,ncols,4)
xx = m.gain_target(~out); yy = y(~out);;
sh = scatter(xx,yy,20,mcv(~out,:));
sh.MarkerFaceColor = 'flat';
sh.MarkerEdgeColor = 'w';
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Slope (PC/dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of slope ~ contrast + gain over mice';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(y,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;


%
% regress threshold against contrast and gain (ADAPTATION)
y = m.thresh;
[~,out] = rmoutliers(y);
%formula = 'thresh ~ contrast + gain_adapt + (1|mouse)';
formula = 'thresh ~ contrast + gain_adapt';
lme = fitlme(m(~out,:),formula);

subplot(nrows,ncols,7)
xx = m.gain_adapt(~out); yy = y(~out);
sh = scatter(xx,yy,20,mcva(~out,:));
sh.MarkerFaceColor = 'flat';
sh.MarkerEdgeColor = 'w';
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Threshold (dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of threshold ~ contrast + gain over mice';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(y,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;


%
% regress slope against contrast and gain (ADAPTATION)
y = m.slope;
[~,out] = rmoutliers(y);
%formula = 'slope ~ contrast + gain_adapt + (1|mouse)';
formula = 'thresh ~ contrast + gain_adapt';
lme = fitlme(m(~out,:),formula);

subplot(nrows,ncols,8)
xx = m.gain_adapt(~out); yy = y(~out);
sh = scatter(xx,yy,20,mcva(~out,:));
sh.MarkerFaceColor = 'flat';
sh.MarkerEdgeColor = 'w';
rh = refline; rh.LineWidth = 1; plotPrefs;
xlabel('Gain'); ylabel('Psychometric Slope (PC/dB)');
title(sprintf('Contrast: %g\nGain: %g',lme.Coefficients.pValue(2:end)))

% stats
statcnt = statcnt + 1;
stat(statcnt).type = 'lme of slope ~ contrast + gain over mice';
stat(statcnt).test = 'fitlme';
stat(statcnt).formula = formula;
stat(statcnt).n = size(y,1);
stat(statcnt).stats = lme;
stat(statcnt).coeffs = lme.Coefficients;