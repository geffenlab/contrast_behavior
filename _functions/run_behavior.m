function res = run_behavior

%addpath(genpath('~/gits/gain-gonogo/_analysis/'));

%%%%%%%%
%% SETUP
%%%%%%%%
% directories/mice to consider
filePath = '~/gits/gain-gonogo/_data';
addpath(genpath('~/gits/gain-gonogo/_analysis/'));
addpath(genpath('./_functions/'));
mouseList = {'CA046','CA047','CA048','CA049','CA051','CA052','CA055',...
             'CA061','CA070','CA072','CA073','CA074','CA075','CA102',...
             'CA104','CA106','CA107',...
             'CA118','CA119','CA121','CA122',...
             'CA123','CA124','CA125','CA126'};

% adjustment for psychometric data
%adjust = '';
%adjust = '_2N';
adjust = '_log';

% results file to use by default
resFile = fullfile('_data',sprintf('_res_behavior_psych_210616%s.mat',adjust));

% compute percent correct (otherwise, does p(respond))
pcFlag = true;

% if results aren't there, make them, otherwise, load them
if ~exist(resFile,'file');
    buildBehaviorRes(mouseList,filePath,resFile,pcFlag,adjust);
end
load(resFile);

% remove dummy mouse
dat(contains(dat.mouse,'CA999'),:) = [];

% important variables
colors = [0 0 1; 1 0 0];
colors_lite = [.7 .7 1; 1 .7 .7];
cc = {'b','r'};
ms = 5;
figSize = [800 1400];
mdl = fitLogGrid();
thresh2use = 'thresh';
groupStatsBy = 'mouse';
%groupStatsBy = 'date';
labels = {'Low Contrast', 'High Contrast'};
axOff = .05;


if strcmp(groupStatsBy,'mouse')
    connectTheDots = true;
else
    connectTheDots = false;
end

% plot name
pname = fullfile('./_plots',...
                 sprintf('_behavior_res_%s_%s.pdf',...
                         thresh2use,...
                         groupStatsBy));

% break out volumes and get the number of sessions/mice per
lgth = cellfun(@length,dat.vols);
tmp = dat(lgth == 6,:);
szz = cellfun(@size,tmp.vols,'uniformoutput',false);
szz = cat(1,szz{:});
tmp = tmp(szz(:,1) == 1,:);
vol6 = cat(1,tmp.vols{:});
[urows,ia,ic] = unique(vol6,'row');

% for each set of volumes, print the volumes, number of mice,
% number of sessions
for i = 1:size(urows,1)
    I = ic == i;
    mice{i,:} = unique(tmp.mouse(I));
    sess(i) = sum(I);
    
    for j = 1:2
        I = ic == i & tmp.contrast == (j-1);
        sess_c(i,j) = sum(I);
    end
end

[~,sortI] = sort(urows(:,end),'descend');

for i = 1:length(sortI)
    ind = sortI(i);
    disp([urows(ind,:), length(mice{ind,:}), sess(ind), sess_c(ind,:)]);
    mm = mice{ind,:}
end

disp(urows(sortI,:))
disp(mice(sortI,:))
disp(sess(sortI))




%% plot training performance for each mouse and mean
% index training sessions and average the performance

fmain = figure(11); clf;
saveFigPDF(fmain,figSize,pname);
nrows = 7;
ncols = 4;

% keep only training sessions, discarding ones with muscimol, or no
% noise condition
include = dat.sesscode == 1 & ~contains(dat.note,'muscimol') & ~contains(dat.label,'NoNoise');
t = dat(include,:);
t.mean_pc = mean(vertcat(t.pc{:}),2);

% plot
ax01 = subplot(nrows,ncols,[1 2]); hold on;
plotTrainingPerf(t,mouseList,colors);
ax01.Position(3) = ax01.Position(3) - axOff;



%%%%%%%%%%%%%%%%%%%%%
%% PSYCHOMETRIC STUFF
%%%%%%%%%%%%%%%%%%%%%
clear t;

% select only psych sessions and clean up cell arrays in the table
include = dat.sesscode == 2 & ~contains(dat.note,'muscimol') & ~contains(dat.label,'NoNoise');
data = dat(include,:);
data = cleanTable(data);

% exclude sessions with high false alarms and very low thresholds
faCut = .5;
threshCut = -inf;
exclude = data.rate(:,1) <= faCut & data.thresh >= threshCut;
data = data(exclude,:);

% temporary table
t = data;

% try using 75% threshold instead
if strcmp(thresh2use,'thresh75');
    t.thresh = t.thresh75;
end


%% mouse means
f1 = figure(10); clf;
[mouseSNR,mousePerf,mouseContrast,mouseID] = plotMouseCurves(t,mouseList,mdl,colors,ms);
saveFigPDF(f1,[800 800],'./_plots/_behavior_psych_curves.pdf',.2)

%% thresholds over time
figure(fmain);
nrows = 7;
ncols = 4;

% get first training date for each mouse
dt0 = varfun(@min,dat,'InputVariables','date',...
             'GroupingVariables','mouse');

% apply to each mouse's psychometric data
t.expDays = zeros(size(t,1),1);
for i = 1:size(dt0,1)
    I = contains(t.mouse,dt0.mouse{i});
    t.expDays(I) = round(days(t.date(I) - dt0.min_date(i)));
    
end

% plot
ax02 = subplot(nrows,ncols,[3]); hold on;
for i = 1:2
    I = t.contrast == (i-1);
    plot(t.expDays(I),t.thresh(I),'.',...
         'Color',colors(i,:)+(colors(i,:)==0)*.6);
    [b,p,xp,yp,ypci] = fitlmWrapper(t.expDays(I),t.thresh(I));
    [h1,h2,h3] = plotlmWrapper(xp,yp,ypci,colors(i,:));
    [psym,pval] = pvalStr(p);
    text(.7,.8+.1*(i-1),sprintf('p=%s %s',pval,psym),...
         'units','normalized','color',colors(i,:))
    
end
xlabel('Time (days rel. task exposure)');
ylabel('Threshold (dB SNR)'); axis tight; plotPrefs;
ax02.Position(2) = ax01.Position(2);
ax02.Position(4) = ax01.Position(4);

%% false alarms
% stats per mouse
mp = grpstats(t,{'mouse','contrast'},'mean',...
              'datavars',{'rate'});
mp = nanpadTable(mp,'mouse','contrast');

% threshold distributions
ax4 = subplot(nrows,ncols,4); hold on;
var = 'mean_rate';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem',...
             'FaceAlpha',.5);
plotPrefs; ylabel('False Alarm Rate');
ylim([0 1]);

% stats
statcnt = 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)(:,1)));
res(statcnt).type = 'paired ttest on fa rate for all mice';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',1);




%%%%%%%%%%%%%%%%%%%%%%%
%% aggregate psych data
%%%%%%%%%%%%%%%%%%%%%%%
xl = [-10 30];
% for each contrast, plot the psychometric curves
figure(fmain); 
ax5 = subplot(nrows,ncols,5); hold on
for i = 1:2
    
    % select data
    I = mouseContrast == (i-1);
    x = round(mouseSNR(I),1);
    y = mousePerf(I);
    m = mouseID(I);
    
    % line per mouse
    um = unique(m);
    for j = 1:length(um)
        mi = m == um(j);
        plot(x(mi),y(mi),'color',[colors(i,:) .15])
    end
        
    
    
    % stats over mice
    mx = grpstats(x,x,'mean');
    my = grpstats(y,x,'mean');
    sy = grpstats(y,x,'sem');
    ny = grpstats(x,x,'numel');
    
    % fit
    xf = linspace(min(mx),max(mx),100);
    %[prms,mdl,thresh,sense,~,~,thresh75] = fitLogGrid(mx,my,[],ny,[],.75);
    [prms,mdl,thresh,sense,~,~,thresh75] = fitLogGrid(x(:),y(:),[],[],[],.75);
    
    if strcmp(thresh2use,'thresh75')
        thresh = thresh75;
    end
    
    % plot
    plot(xl,[.5 .5],'k');
    %errorBars(mx',my',colors_lite(i,:),'sem',[],sy,[],...
%          'o','LineWidth',1,'MarkerSize',ms, ...
%           'MarkerFaceColor',colors_lite(i,:));
% errorbar(mx,my,sy,'o','MarkerSize',ms,'CapSize',0,...
%         'Color',colors(i,:)+(colors(i,:)==0)*.7,...
%         'FaceColor',colors(i,:)+(colors(i,:)==0)*.7);
    plot(xf,mdl(prms,xf),'Color',colors(i,:),'LineWidth',1);
    plot([thresh thresh],[.5 mdl(prms,thresh)],'--','Color',colors(i,:))
    
end  
xlim(xl);
ylabel('Percent Correct'); plotPrefs;


% stats per mouse
mp = grpstats(t,{'mouse','contrast'},'mean',...
              'datavars',{'thresh','maxslope','sense'});
mp = nanpadTable(mp,'mouse','contrast');

% threshold distributions
ax6 = subplot(nrows,ncols,6); hold on;
var = 'mean_thresh';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem',...
             'FaceAlpha',.5);
plotPrefs; ylabel('Threshold (dB SNR)');


% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on thresholds for all mice';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',1);



% slope
ax7 = subplot(nrows,ncols,7); hold on;
var = 'mean_maxslope';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem',...
             'FaceAlpha',.5);
plotPrefs; ylabel('Slope (PC/dB)');

% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on slopes for all mice';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',1);

% sensitivity
ax8 = subplot(nrows,ncols,8); hold on;
var = 'mean_sense';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem',...
             'FaceAlpha',.5);
plotPrefs; ylabel('Sensitivity (\beta)');

% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on sensitivity for all mice';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',1);

drawnow;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comparison purely low and high contrast %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set up handles for different conditions
t.rangeI = range(t.vols,2);
t.meanI = mean(t.vols,2);
t.pequalI = all(diff(t.volP(:,2:end),1,2)==0,2);
t.requalI = all(diff(t.rewcont(:,2:end),1,2)==0,2);

% look at unique combinations of them
mInd = [t.contrast t.meanI t.rangeI t.pequalI t.requalI];
uSess = unique(mInd,'row');


%% different contrasts, same mean and range
% session selection
selector = t.meanI == 12.5 & t.rangeI == 25 & ...
    t.pequalI == 1 & t.requalI == 1;

selector = t.rangeI == 25 & t.requalI == 1;

labels = {'Low Contrast', 'High Contrast'};

% stats (exclude mice with fewer than 3 sessions)
mp = grpstats(t(selector,:),{groupStatsBy,'contrast'},'mean',...
              'datavars',{'thresh','maxslope','vols','pc','sense'});
mp = nanpadTable(mp,groupStatsBy,'contrast',3);

% exclude 

% curve plot
ax9 = subplot(nrows,ncols,[9]); hold on;
plot(xl,[.5 .5],'k');
for i = 1:2
    I = mp.contrast == (i-1);
    plotPsychCurves(mp,I,colors(i,:),'-',colors_lite(i,:),ms,thresh2use);
end
ylabel('Percent Correct'); 
plotPrefs; axis tight; ylim([.4 1]); xlim(xl);
  

% threshold
ax10 = subplot(nrows,ncols,10); hold on;
var = 'mean_thresh';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem','FaceAlpha',.5);
plotPrefs; ylabel('Threshold (dB SNR)');

% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on thresh for matched targets';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);
res(statcnt).xx = {x1, x2};


% slope
ax11 = subplot(nrows,ncols,11); hold on;
var = 'mean_maxslope';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem','FaceAlpha',.5);
plotPrefs; ylabel('Slope (PC/dB)');

% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on slope for matched targets';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);
res(statcnt).xx = {x1, x2};



% sensitivity
ax12 = subplot(nrows,ncols,12); hold on;
var = 'mean_sense';
x1 = mp.(var)(mp.contrast==0); 
x2 = mp.(var)(mp.contrast==1);
[ph,pl,pe] = plotDist([1 2],[x1 x2],'k',true,[],[],[],'MarkerSize',ms);
barWithError([x1 x2],{'Low Contrast','High Contrast'},{'b','r'},'sem','FaceAlpha',.5);
plotPrefs; ylabel('Sensitivity (\beta)');

% stats
statcnt = statcnt + 1;
[~,pv,~,stats] = ttest(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'paired ttest on sensitivity for matched targets';
res(statcnt).test = 'ttest';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = mean([x1,x2],1,'omitnan');
res(statcnt).std = std([x1,x2],1,'omitnan');
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);
res(statcnt).xx = {x1, x2};


drawnow;


%% low contrast, different ranges, same mean
% session selection
ranges = {15,25};
cselect = 1;
selector = t.contrast == cselect-1 & t.meanI == 12.5 & ...
    t.requalI == 1;

% plot stuff
mfill = {'none',colors_lite(cselect,:)};
ls = {'--','-'}; labels = {'Narrow','Wide'};
alphas = [.1 .5];

% stats
clear mp st;
st = t(selector & (t.rangeI == 15 | t.rangeI==25),:);
grpVars = {groupStatsBy,'rangeI'};
mp = grpstats(st,grpVars,'mean',...
              'datavars',{'thresh','maxslope','vols','pc','sense'});


% curve plots
ax13 = subplot(nrows,ncols,[13]); hold on;
plot(xl,[.5 .5],'k');
for i = 1:2
    I = mp.rangeI == ranges{i};
    plotPsychCurves(mp,I,colors(cselect,:),ls{i},mfill{i},ms,thresh2use,false);
end
ylabel('Percent Correct'); 
plotPrefs; axis tight; ylim([.4 1]); xlim(xl);

    
% threshold distributions
ax14 = subplot(nrows,4,14); hold on;
var = 'mean_thresh'; grp = 'rangeI'; grp_val = ranges;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Threshold (dB SNR)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val{1}); 
x2 = mp.(var)(mp.(grp)==grp_val{2});
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on thresh low contrast range=[15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);



% slope
ax15 = subplot(nrows,ncols,15); hold on;
var = 'mean_maxslope'; grp = 'rangeI'; grp_val = ranges;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Slope (PC/dB)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val{1}); 
x2 = mp.(var)(mp.(grp)==grp_val{2});
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on slope low contrast range=[15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);

% sense
ax16 = subplot(nrows,ncols,16); hold on;
var = 'mean_sense'; grp = 'rangeI'; grp_val = ranges;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Sensitivity (\beta)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val{1}); 
x2 = mp.(var)(mp.(grp)==grp_val{2});
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on sense low contrast range=[15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);


drawnow;



%% high contrast, different ranges, same mean
% session selection
ranges = {12,15,25};
cselect = 2;
selector = t.contrast == cselect-1 & ...
    t.requalI == 1;

% plot stuff
mfill = {'none',colors_lite(cselect,:)};
ls = {'--','-'}; labels = {'Narrow','Wide'};
alphas = [.1 .5];

% stats
clear mp st;
st = t(selector & (t.rangeI == 12 | t.rangeI == 15 | t.rangeI==25),:);
st.rangeCmp = st.rangeI == 25 & st.rangeI ~= 15; % test 12/15 vs 25
grpVars = {groupStatsBy,'rangeCmp'};
mp = grpstats(st,grpVars,'mean',...
              'datavars',{'thresh','maxslope','vols','pc','sense','rangeI'});


% psych plot
ax17 = subplot(nrows,4,[17]); hold on;
plot(xl,[.5 .5],'k');
ur = unique(mp.rangeCmp);
for i = 1:length(ur)
    I = mp.rangeCmp == ur(i);
    plotPsychCurves(mp,I,colors(cselect,:),ls{i},mfill{i},ms,thresh2use,false);
end
ylabel('Percent Correct'); 
plotPrefs; axis tight; ylim([.4 1]); xlim(xl);
    
% threshold distributions
ax18 = subplot(nrows,ncols,18); hold on;
var = 'mean_thresh'; grp = 'rangeCmp'; grp_val = ur;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Threshold (dB SNR)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val(1)); 
x2 = mp.(var)(mp.(grp)==grp_val(2));
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on thresh high contrast range=[12|15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);

% slope
ax19 = subplot(nrows,ncols,19); hold on;
var = 'mean_maxslope'; grp = 'rangeCmp'; grp_val = ur;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Slope (PC/dB)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val(1)); 
x2 = mp.(var)(mp.(grp)==grp_val(2));
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on slope high contrast range=[12|15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);

% sense
ax20 = subplot(nrows,ncols,20); hold on;
var = 'mean_sense'; grp = 'rangeCmp'; grp_val = ur;
plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
               alphas,ls,'k',ms,labels,connectTheDots);
plotPrefs; ylabel('Sensitivity (\beta)');

% stats
statcnt = statcnt + 1;
x1 = mp.(var)(mp.(grp)==grp_val(1)); 
x2 = mp.(var)(mp.(grp)==grp_val(2));
[~,pv,~,stats] = ttest2(x1,x2);
plotPvalSmart(pv,[1 2],max(mp.(var)));
res(statcnt).type = 'unpaired ttest on sense high contrast range=[12|15,25]';
res(statcnt).test = 'ttest2';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).means = [mean(x1,'omitnan'),mean(x2,'omitnan')];
res(statcnt).std = [std(x1,'omitnan'),std(x2,'omitnan')];
res(statcnt).n = [sum(~isnan(x1)) sum(~isnan(x2))];
res(statcnt).effectSizeMethod = 'Hedges g';
res(statcnt).effectSize = mes(x1,x2,'hedgesg','isDep',0);

drawnow;





%%%%%%%%%%%%%%%%%%
%% OFFSET STUFF %%
%%%%%%%%%%%%%%%%%%

offResFile = fullfile('_data','_res_behavior_offset_210428.mat');
[t,uOff] = buildOffsetBehaviorRes(dat,offResFile);


%% offset analysis
% include only sessions with low FA rates
include = mean(t.rate_0,2,'omitnan') < .5;

gm = grpstats(t(include,:),{'mouse','contrast'},'mean','DataVars',{'offs_0','pc_0','tau'});
xf = linspace(min(uOff),max(uOff),100);

% combination index for adjacent times
clear ind;
ind{1} = [1 2];
ind{2} = [3 4];
ind{3} = [5 6];
ind{4} = [7 8];
ind{5} = [9 10];
offs = combindstats(mean(gm.mean_offs_0,'omitnan'),ind);

% for each mouse, get the mean performance per contrast and fit
% exponential function (at ***threshold***)
uM = unique(gm.mouse);
dm = nan(length(uM),length(ind),2);
taus = nan(length(uM),2);
lags = 1:4;
for i = 1:length(uM)
    for j = 1:2
        I = contains(gm.mouse,uM{i}) & gm.contrast == (j-1);
        if sum(I) > 0
            x = offs; y = combindstats(gm.mean_pc_0(I,:),ind);
            dm(i,:,j) = y;
            [params,mdl,taus(i,j)] = fitExpGrid(x(lags),y(lags));
        end
    end
end

% check for significant changes in discriminability per contrast
% and timepoint
statcnt = statcnt + 1;
clear stats;
for j = 1:2
    for i = 1:size(dm,2)-1
        [pv(j,i),~,stats{j,i}] = signrank(squeeze(dm(:,1,j)), ...
                                         squeeze(dm(:,i+1,j)));
        effsz(j,i) = stats{j,i}.zval / sqrt(size(dm,1));
    end
end
res(statcnt).type = 'signrank discriminability: 1st time against others';
res(statcnt).test = 'signrank';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).medians = squeeze(median(dm,1,'omitnan'))';
res(statcnt).iqr = squeeze(iqr(dm,1))';
res(statcnt).n = squeeze(sum(~isnan(dm),1))';
res(statcnt).effectSizeMethod = 'Z/sqrt(n)';
res(statcnt).effectSize = effsz;
res(statcnt).h_fdr = fdr_bh(pv);
        
% plot
ax33 = subplot(nrows,ncols,nrows*ncols-[2 3]); hold on;
for i = 1:2
    x = offs; y = squeeze(dm(:,:,i));;
    xf = linspace(min(x),max(x),100);
    
    % fit
    xx = x(~isnan(x(:))); yy = y(~isnan(y(:)));
    [params,mdl,tau(i)] = fitExpGrid(x,mean(y,'omitnan'));
    plot(xf,mdl(params,xf),'color',colors(i,:),'linewidth',1);
    
    
    % plot
    e = errorBars(x,y,colors_lite(i,:));
    e.Marker = 'o'; e.LineStyle = 'none';
    e.MarkerFaceColor = colors_lite(i,:);
    e.MarkerSize = ms;
    

end
YL = ylim;

% plot test results
cnt = 0;
for k = 1:2
    for ii = 1:4
        if res(statcnt).h_fdr(k,ii)
            cnt = cnt + 1;
            plot([x(1) x(ii+1)],...
                 [YL(2)+.05*cnt*diff(YL) YL(2)+.05*cnt*diff(YL)],...
                 'color',colors(k,:),'linewidth',.75)
        end
    end
end
plot([0 0],ylim,'k--');
axis tight; xlim([-.1 1.1]);
xlabel('Time (s)');
ylabel('Percent Correct');
plotPrefs;


%% taus
ax35 = subplot(nrows,ncols,nrows*ncols - 1); hold on;
[ph,pl,pe] = plotDist([1 2],taus,'k',true,[],[],[],'MarkerSize',ms);
barWithError(taus,{'Low Contrast','High Contrast'},{'b','r'},'sem','FaceAlpha',.5);
ylabel('\tau (s)'); plotPrefs;

% stats
statcnt  = statcnt + 1;
[pv,hh,stats] = ranksum(taus(:,1),taus(:,2));
plotPvalSmart(pv,[1 2],max(taus(:)));
res(statcnt).type = 'ranksum test of taus between low and high';
res(statcnt).test = 'ranksum';
res(statcnt).p = pv;
res(statcnt).stats = stats;
res(statcnt).medians = median(taus,'omitnan');
res(statcnt).iqr = iqr(taus);
res(statcnt).n = sum(~isnan(taus));
res(statcnt).effectSizeMethod = 'Z/sqrt(n)';
res(statcnt).effectSize = stats.zval / sqrt(size(taus,1));
res(statcnt).xx = taus;


if false
    % old version of analysis
    gm = grpstats(t,{'mouse','contrast'},'mean','DataVars',{'offs_0','pc_0','tau'});
    xf = linspace(min(uOff),max(uOff),100);

    % fit the mean for each mouse
    clear mtau;
    for i = 1:size(gm,1)
        
        x = gm.mean_offs_0(i,:); y = gm.mean_pc_0(i,:);
        [params,mdl,mtau(i)] = fitExpGrid(x,y);
        
    end
    gm.mtau = mtau';

    % save copy of all mice data
    allm = gm;

    % include only mice with >1 session and no early targets
    include = gm.GroupCount>1 & ~isnan(gm.mean_offs_0(:,2));
    gm = gm(include,:);

    ax33 = subplot(nrows,ncols,nrows*ncols - [2 3]); hold on;
    plot([0 0],[.6 .9],'k--');
    for i = 1:2
        I = gm.contrast == (i-1);
        x = gm.mean_offs_0(I,:); y = gm.mean_pc_0(I,:);
        
        % plot
        e = errorBars(nanmean(x),y,colors_lite(i,:));
        e.Marker = 'o'; e.LineStyle = 'none';
        e.MarkerFaceColor = colors_lite(i,:);
        e.MarkerSize = ms;
        
        % fit
        xx = x(~isnan(x(:))); yy = y(~isnan(y(:)));
        [params,mdl,tau(i)] = fitExpGrid(xx,yy);
        plot(xf,mdl(params,xf),'color',colors(i,:),'linewidth',1);
    end
    xlim([-.1 1.1]); ylim([.6 .9])
    xlabel('Time (s., rel. contrast switch)');
    ylabel('Percent Correct');
    plotPrefs;


    ax35 = subplot(nrows,ncols,nrows*ncols - 1); hold on;
    var = 'mtau'; grp = 'contrast'; grp_val = {0, 1};
    labels = {'Low Contrast','High Contrast'};
    connectTheDots = true;
    [h,b] = plotPsychStats(gm,var,grp,grp_val,colors,...
                           [.5 .5],{'-','-'},'k',ms,labels,connectTheDots);
    pv = ranksum(mtau(gm.contrast==0),mtau(gm.contrast==1));
    plotPvalSmart(pv,[1 2],max(gm.(var)));
    ylabel('\tau (s)'); plotPrefs;
end


%% shift all axes inwards
ax4.Position(3) = ax4.Position(3) - .05;

ax6.Position(3) = ax6.Position(3) - .05;
ax7.Position(3) = ax7.Position(3) - .05;
ax8.Position(3) = ax8.Position(3) - .05;

ax10.Position(3) = ax10.Position(3) - .05;
ax11.Position(3) = ax11.Position(3) - .05;
ax12.Position(3) = ax12.Position(3) - .05;

ax14.Position(3) = ax14.Position(3) - .05;
ax15.Position(3) = ax15.Position(3) - .05;
ax16.Position(3) = ax16.Position(3) - .05;

ax18.Position(3) = ax18.Position(3) - .05;
ax19.Position(3) = ax19.Position(3) - .05;
ax20.Position(3) = ax20.Position(3) - .05;

ax33.Position(3) = ax33.Position(3) - .05;
ax35.Position(3) = ax35.Position(3) - .05;



%% clean up axes
ax4.Position([2,4]) = ax02.Position([2,4]);

ax6.Position([2,4]) = ax5.Position([2,4]);
ax7.Position([2,4]) = ax5.Position([2,4]);
ax8.Position([2,4]) = ax5.Position([2,4]);

ax10.Position([2,4]) = ax9.Position([2,4]);
ax11.Position([2,4]) = ax9.Position([2,4]);
ax12.Position([2,4]) = ax9.Position([2,4]);

ax14.Position([2,4]) = ax13.Position([2,4]);
ax15.Position([2,4]) = ax13.Position([2,4]);
ax16.Position([2,4]) = ax13.Position([2,4]);

ax18.Position([2,4]) = ax17.Position([2,4]);
ax19.Position([2,4]) = ax17.Position([2,4]);
ax20.Position([2,4]) = ax17.Position([2,4]);

ax35.Position(2) = ax33.Position(2);
ax35.Position(4) = ax17.Position(4);

saveFigPDF(fmain,figSize,pname,.2);

    
































if (2 + 2) == 5
    %% low contrast, diff mean, same range
    cselect = 1;
    selector = t.contrast == cselect-1 & t.rangeI == 25;

    % plot stuff
    mfill = {'none',colors_lite(cselect,:)};
    ls = {'--','-'}; labels = cellstr(num2str([means{:}]'))';
    alphas = [.1 .5];

    % stats
    clear mp st;
    st = t(selector & (t.meanI == means{1} | t.meanI==means{2}),:);
    mp = grpstats(st,{groupStatsBy,'meanI'},'mean',...
                  'datavars',{'thresh','maxslope','vols','pc'});

    % psych plot
    ax21 = subplot(9,4,[21]); hold on;
    plot([-5 30],[.5 .5],'k');
    for i = 1:length(means)
        I = mp.meanI == means{i};
        plotPsychCurves(mp,I,colors(cselect,:),ls{i},mfill{i},ms,thresh2use);
    end
    ylabel('Percent Correct'); 
    plotPrefs; axis tight; ylim([.4 1]); xlim([-5 30]);



    % threshold distributions
    ax22 = subplot(9,4,22); hold on;
    var = 'mean_thresh'; grp = 'meanI'; grp_val = means;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Threshold (dB SNR)');



    % slope
    ax23 = subplot(9,4,23); hold on;
    var = 'mean_maxslope'; grp = 'meanI'; grp_val = means;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Slope (PC/dB)');





    %% high contrast, reward contingency
    rewC = {0 1};
    cselect = 2;
    selector = t.contrast == cselect-1 & t.rangeI == 25 & ...
        t.meanI == 12.5 & t.pequalI == 0;

    % plot stuff
    mfill = {'none',colors_lite(cselect,:)};
    ls = {'--','-'}; labels = {'Biased Rewards', 'Equal Rewards'};
    alphas = [.1 .5];

    % stats
    clear mp st;
    st = t(selector & (t.requalI==rewC{1} | t.requalI==rewC{2}),:);
    mp = grpstats(st,{groupStatsBy,'requalI'},'mean',...
                  'datavars',{'thresh','maxslope','vols','pc'});

    % psych plot
    ax25 = subplot(9,4,[25]); hold on;
    plot([-5 30],[.5 .5],'k');
    for i = 1:length(means)
        I = mp.requalI == rewC{i};
        plotPsychCurves(mp,I,colors(cselect,:),ls{i},mfill{i},ms,thresh2use);
    end
    xlabel('Target Level (dB SNR)'); ylabel('Percent Correct'); 
    plotPrefs; axis tight; ylim([.4 1]); xlim([-5 30]);



    % threshold distributions
    ax26 = subplot(9,4,26); hold on;
    var = 'mean_thresh'; grp = 'requalI'; grp_val = rewC;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Threshold (dB SNR)');



    % slope
    ax27 = subplot(9,4,27); hold on;
    var = 'mean_maxslope'; grp = 'requalI'; grp_val = rewC;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Slope (PC/dB)');




    %% high contrast, presentation probability
    presP = {0 1};
    cselect = 2;
    selector = t.contrast == cselect-1 & t.rangeI == 25 & ...
        t.meanI == 12.5 & t.requalI == 1;

    % plot stuff
    mfill = {'none',colors_lite(cselect,:)};
    ls = {'--','-'}; labels = {'Biased Targets', 'Equal Targets'};
    alphas = [.1 .5];

    % stats
    clear mp st;
    st = t(selector & (t.pequalI == presP{1} | t.pequalI==presP{2}),:);
    mp = grpstats(st,{groupStatsBy,'pequalI'},'mean',...
                  'datavars',{'thresh','maxslope','vols','pc'});

    % psych plot
    ax29 = subplot(9,4,[29]); hold on;
    plot([-5 30],[.5 .5],'k');
    for i = 1:length(means)
        I = mp.pequalI == presP{i};
        plotPsychCurves(mp,I,colors(cselect,:),ls{i},mfill{i},ms,thresh2use);
    end
    xlabel('Target Level (dB SNR)'); ylabel('Percent Correct'); 
    plotPrefs; axis tight; ylim([.4 1]); xlim([-5 30]);

    % threshold distributions
    ax30 = subplot(9,4,30); hold on;
    var = 'mean_thresh'; grp = 'pequalI'; grp_val = presP;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Threshold (dB SNR)');


    % slope
    ax31 = subplot(9,4,31); hold on;
    var = 'mean_maxslope'; grp = 'pequalI'; grp_val = presP;
    plotPsychStats(mp,var,grp,grp_val,colors(cselect,:),...
                   alphas,ls,'k',ms,labels,connectTheDots)
    [~,pv] = ttest2(mp.(var)(mp.(grp)==grp_val{1}), ...
                    mp.(var)(mp.(grp)==grp_val{2}));
    plotPvalSmart(pv,[1 2],max(mp.(var)));
    plotPrefs; ylabel('Slope (PC/dB)');
    
    ax22.Position(3) = ax22.Position(3) - .05;
    ax23.Position(3) = ax23.Position(3) - .05;
    ax26.Position(3) = ax26.Position(3) - .05;
    ax27.Position(3) = ax27.Position(3) - .05;
    ax30.Position(3) = ax30.Position(3) - .05;
    ax31.Position(3) = ax31.Position(3) - .05;
    
    
    ax22.Position(2) = ax21.Position(2);
    ax22.Position(4) = ax21.Position(4);
    ax23.Position(2) = ax21.Position(2);
    ax23.Position(4) = ax21.Position(4);
    ax26.Position(2) = ax25.Position(2);
    ax26.Position(4) = ax25.Position(4);
    ax27.Position(2) = ax25.Position(2);
    ax27.Position(4) = ax25.Position(4);
    ax30.Position(2) = ax29.Position(2);
    ax30.Position(4) = ax29.Position(4);
    ax31.Position(2) = ax29.Position(2);
    ax31.Position(4) = ax29.Position(4);
end




        
        
        





    
    

