function run_normative_model_inactivation(sigma_n,sigma_c,gc,f,silence)

close all;

%% initialize 
fprintf('--------------\nRunning normative model simulations...\n');

% inputs
if nargin < 1
    % neural noise level (1 = unit sigma)
    sigma_n = 1;
end

if nargin < 2
    % contrast noise levels
    sigma_c = [1 2];
end

if nargin < 3
    % gain control ('off' means average gain, 'low' means low gain,
    % [] means full gain control)
    gc = [];
end

if nargin < 4
    % target noise scaling factor (sigma_target = sigma_c / f)
    f = .25;
end

if nargin < 5
    % background scaling factor (sigma_background = sigma_background*silence)
    silence = 1;
end

%define means of target distribution
muTlow = 0:.25:3;

% contrast widths
sigmaLow = sigma_c(1);
sigmaHigh = sigma_c(2);

%index of distribution used in previous plots
indTarget = find(muTlow == 1.5);


%% generate results for all values of muT
tag = sprintf('_sigman-%d_sigmac-%g-%g_gc%s_f-%g_silence-%g',...
              sigma_n,sigma_c(1),sigma_c(2),gc,f,silence);
resFile = fullfile('./_data',['normative_model' tag '.mat']);

% plot location
plot_dir = sprintf('./_plots/%s',tag);
mkdir(plot_dir)

%% plot figures used in Cosyne abstract
[res,p] = varianceSims(muTlow(indTarget),sigmaLow,sigmaHigh,f,sigma_n,gc,silence);
plotFigs_CFA(res,plot_dir); drawnow;



%% generate res for all target means
if ~exist(resFile)
    fprintf('\nSimulating different target volumes...\n');
    for i = 1:length(muTlow)
        rr(i) = varianceSims(muTlow(i),sigmaLow,sigmaHigh,f,sigma_n,gc,silence);
    end
    save(resFile,'rr','-v7.3');
else
    fprintf('Loading simulation results...\n');
    load(resFile);
end




%% compare performance metrics
nLevels = 16;
times = 1:50;
iT = 1:length(times);
bins = -.5:(nLevels-.5);
bb = [bins,nLevels+.5,nLevels+1.5]-1;

% precompute distributions for each contrast, volume, time
for j = 1:length(rr)
    for i=1:numel(iT)
        
        % distributions
        targetDist(:,j,i,1) = rr(j).buffer.yT(50+times(iT(i)),:);
        targetDist(:,j,i,2) = rr(j).buffer.yT(100+times(iT(i)),:);
        noiseDist(:,j,i,1) = rr(j).buffer.yB(50+times(iT(i)),:);
        noiseDist(:,j,i,2) = rr(j).buffer.yB(100+times(iT(i)),:);
        
        % during each run, what was the gain?
        gain(:,j,i,1) = rr(j).buffer.k(50+times(iT(i)),:);
        gain(:,j,i,2) = rr(j).buffer.k(100+times(iT(i)),:);
        
    end
end

% decision rule: criterion that optimizes performance
% across all volumes and times
for i = 1:2
    allT = targetDist(:,:,:,i);
    allB = noiseDist(:,:,:,i);
    
    allVals = [allT(:); allB(:)];
    trueVals = [ones(size(allT(:))); zeros(size(allB(:)))];
    
    % try criterion to optimize detection
    crits = 0:16;
    for j = 1:length(crits)
        pred = allVals > crits(j);
        pc(j) = mean(pred == trueVals);
    end
    [~,mxi] = max(pc);
    crit(i) = crits(mxi);
    
end

contrastTimes = [50 100];
for j = 1:length(rr)
    for i=1:numel(iT)
        for k = 1:2
            
            % distributions for this condition
            tDist = squeeze(targetDist(:,j,i,k));
            nDist = squeeze(noiseDist(:,j,i,k));

            % roc values
            [~,~,roc(j,i,k)] = computeROC(nDist,tDist);
            
            % KL divergence
            dkl(j,i,k) = rr(j).buffer.Dkl(contrastTimes(k)+times(iT(i)),2);
            
            % distribution overlap
            doverlap(j,i,k) = rr(j).buffer.Dkl(contrastTimes(k)+times(iT(i)),1);
            
            % criterion decision rule
            allResp = [tDist; nDist];
            labels = [ones(size(tDist)); zeros(size(nDist))];
            pred = allResp > crit(1);
            fa(j,i,k) = mean(pred(labels==0));
            hr(j,i,k) = mean(pred(labels==1));
            perf(j,i,k) = mean(pred == labels);
            
        end
    end
end


% fit curves to each psychometric function
fprintf('Fitting model psychometric functions... '); tic;
metrics = {'roc','dkl','doverlap','perf'};
for m = 1:length(metrics)
    for i = 1:length(iT)
        for j = 1:2
            
            x = muTlow';
            eval(sprintf('y = squeeze(%s(:,i,j));',metrics{m}));
            [prms(m,i,j,:), mdl, thresh(m,i,j), sense(m,i,j), ~, ~,thresh75(m,i,j)] = ...
                fitLogGrid(x,y,[],[],[],.75);
            slope(m,i,j) = max(diff(y)./diff(x));
            
        end
    end
end
toc;

            


% plot graphs for each metric
tI = 1:25;
tcolors{1} = [ones(length(tI),1) ...
              linspace(.8,0,length(tI))' ...
              linspace(.8,0,length(tI))'];
tcolors{2} = [linspace(.8,0,length(tI))'...
              linspace(.8,0,length(tI))'...
              ones(length(tI),1)];
mmarks = linspace(.25,3,length(muTlow));
g = squeeze(mean(gain,[1 2]));


% colors
cols = [1 0 0; 0 0 1];
cols_lite = [1 .6 .6;
             .6 .6 1];

%plot discriminability versus target mean
f2 = figure(21); clf; hold on; 
mi = find(contains(metrics,'doverlap'));
times = 1:5:26;
for i=1:numel(times)
    subplot(3,numel(times),i);hold on;
    x = muTlow'; xf = linspace(min(x),max(x),100);
    
    % high contrast
    y = doverlap(:,times(i),1);
    pp(i,1,:) = squeeze(prms(mi,times(i),1,:));
    tr(i,1) = thresh(mi,times(i),1);
    sl(i,1) = slope(mi,times(i),1);
    plot(x,y,'.','color',cols_lite(1,:),'markersize',10);
    plot(xf,mdl(pp(i,1,:),xf),'color',cols(1,:),'linewidth',1);
    plot([tr(i,1) tr(i,1)],[.5 mdl(pp(i,1,:),tr(i,1))],...
         '--','color',cols(1,:),'linewidth',1);
    
    % low contrast
    y = doverlap(:,times(i),2);
    pp(i,2,:) = squeeze(prms(mi,times(i),2,:));
    tr(i,2) = thresh(mi,times(i),2);
    sl(i,2) = slope(mi,times(i),2);
    plot(x,y,'.','color',cols_lite(2,:),'markersize',10);
    plot(xf,mdl(pp(i,2,:),xf),'color',cols(2,:),'linewidth',1);
    plot([tr(i,2) tr(i,2)],[.5 mdl(pp(i,2,:),tr(i,2))],...
         '--','color',cols(2,:),'linewidth',1);
    xlabel('target mean');
    ylabel('discriminability');
    ylim([.5,1]);xlim([min(muTlow),max(muTlow)])
    plotPrefs;
    
    subplot(3,numel(times),numel(times)+i);
    b = bar([2 1],tr(i,:),'FaceColor','flat');
    b.CData(1,:) = cols(2,:);
    b.CData(2,:) = cols(1,:);
    ylabel('Threshold'); ylim([0 2.5]);
    plotPrefs
    
    subplot(3,numel(times),numel(times)*2+i);
    b = bar([2 1],sl(i,:),'FaceColor','flat');
    b.CData(1,:) = cols(2,:);
    b.CData(2,:) = cols(1,:);
    ylabel('Slope'); ylim([0 .2]);
    plotPrefs;
    
    
end

saveFigPDF(f2,[900 400],fullfile(plot_dir,'_psychCurvesDT.pdf'));


%% subplot for each target mean over time
f11 = figure(11); clf; hold on;

% color map
cw = [1,1,1];
cr = [1,0,0];
cb = [0,0,1];
cRL = [];cBL = [];
for i=1:15
    cRL = [cRL;cw+(cr-cw)*i./(numel(muTlow)+3)];
    cBL = [cBL;cw+(cb-cw)*i./(numel(muTlow)+3)];
end
cmap(:,:,1) = cBL;
cmap(:,:,2) = cRL;

for i = 1:numel(muTlow)
    
    subplot(4,4,i); hold on;
    clear d;
    d(1,:) = [rr(i).buffer.Dkl(95:100,1);rr(i).buffer.Dkl(1:50,1)];
    d(2,:) = [rr(i).buffer.Dkl(45:50,1);rr(i).buffer.Dkl(51:100,1)];
   
    plot([0 0],[.5 .9],'k--');
    title(sprintf('Target mu = %3.2f',muTlow(i)));
    
    for j = 1:size(d,1)
        
        % plot
        clear DD;
        DD = d(j,:);
        plot(-5:50,DD,'linewidth',2,'color',cmap(i+2,:,j));
        ylim([.5 1]);
        
    end
    plotPrefs;
        
end
saveFigPDF(f11,[800 700],fullfile(plot_dir,'_discriminabilityDTperMu.pdf'));






%% nonlinearity slopes over time
times = 1:50;
nT = numel(times);
x0 = linspace(-10,10,100);
ind = 4;
sigmaNoise = p.rs_noise01.sigmaNoise; % noise added to spike
                                      % response
nLevels = 16; % max spikes per time bin?
tI = [100 50]; % contrast index over time

% make nonlinearity at each timepoint for plotting
clear rV;
for i=1:nT
    
    t = times(i);
    
    % low variance
    [~,rV(1,i,:)] = encodeWithDiscrete(x0, mean(res.buffer.k(100+t,:),2), ...
                                       mean(res.buffer.x0(100+t,:),2), ...
                                       0,nLevels);
    % high variance
    [~,rV(2,i,:)] = encodeWithDiscrete(x0, mean(res.buffer.k(100+t,:),2), ...
                                       mean(res.buffer.x0(100+t,:),2), ...
                                       0,nLevels);
    
end


%% gain plots (k)

f12 = figure(12); clf;
cc = {'b','r'};
ksamps = 40;
for i = 1:2
    
    subplot(2,3,2:3); hold on;
    ts = mean(res.buffer.k(tI(i):tI(i)+ksamps,:),2);
    
    %if i == 2
%    ts = -ts;
        % end
    %ts = normalize(ts,'range');
    
    plot(ts,'color',cc{i},'LineWidth',1);

    
    % find half 'max'
    halfMax = range(ts)/2 + min(ts);
    [~,mi] = min(abs(ts - halfMax));
    plot([mi mi],[0 ts(mi)],'--','color',cc{i},'linewidth',1);
    plot(mi,ts(mi),'o','color',cc{i},'linewidth',1);
    xlabel('Time (samples)');
    ylabel('Gain'); 
    plotPrefs; axis tight; ylim([.2 1]);
    
    halfRange(i) = mi;
    
end

subplot(2,3,6);
b = bar([1 2],halfRange,'FaceColor','flat');
b.CData(1,:) = cols(2,:);
b.CData(2,:) = cols(1,:);
ylabel('Timesteps to Half-Range'); ylim([0 11]);


% plot nonlinearities
subplot(2,3,1); hold on;
mRV = squeeze(rV(1,:,:));
times = 1:15;
nt = length(times);
alphas = linspace(.1,1,nt);
for i = 1:nt
    ph = plot(x0,mRV(i,:),'b');
    ph.Color(4) = alphas(i);
end
xlabel('signal'); ylabel('response');
plotPrefs;


% plot nonlinearities
subplot(2,3,4); hold on;
mRV = squeeze(rV(2,:,:));
for i = 1:nt
    ph = plot(x0,mRV(i,:),'r');
    ph.Color(4) = alphas(i);
end
xlabel('signal'); ylabel('response');
plotPrefs;
saveFigPDF(f12,[450 275],fullfile(plot_dir,'deltaGain.pdf'));



%% summary
mu_lo = 1.5;
mu_hi = 2.25;
cols = [0 0 1; 1 0 0];
cols_lite = [.8 .8 1; 1 .8 .8];

% extract each time course
dd(:,1) = [rr(muTlow == mu_lo).buffer.Dkl(95:100,1);...
          rr(muTlow == mu_lo).buffer.Dkl(1:50,1)];
dd(:,2) = [rr(muTlow == mu_hi).buffer.Dkl(45:50,1);...
          rr(muTlow == mu_hi).buffer.Dkl(51:100,1)];
kk(:,1) = mean([rr(muTlow == mu_lo).buffer.k(95:100,:);...
                rr(muTlow == mu_lo).buffer.k(1:50,:)],2);
kk(:,2) = mean([rr(muTlow == mu_hi).buffer.k(45:50,:);...
                rr(muTlow == mu_hi).buffer.k(51:100,:)],2);

% time sampling, etc
t = -5:50;
tI = t>=0;
p0 = [.7 .1 4];
yy = dd;

% plot discriminability and fit
f12 = figure(12); clf;
subplot(3,3,[1 2]); hold on
for i = 1:2
    [prm(:,i),emdl,tau(i)] = fitExpGrid(t(tI),movmean(yy(tI,i),2),p0);
    plot(t,yy(:,i),'.','markersize',10,'color',cols_lite(i,:));
    plot(t(tI),emdl(prm(:,i),t(tI)),'color',cols(i,:),'linewidth',1);
end
axis tight; ylim([.5 1]);
plot([0 0],ylim,'k--'); plotPrefs;
ylabel('Discriminability'); xlabel('Time (steps)');

% plot gain
subplot(3,3,[4 5]); hold on
for i = 1:2
    plot(t,kk(:,i),'-','color',cols(i,:));
end
axis tight;
plot([0 0],ylim,'k--');
ylabel('Gain');
xlabel('Time (steps)');
plotPrefs;

% plot psychometric function near full adaptation
clear sl pp tr;
adaptT = 26;
mi = 3;
mdl = fitLogGrid();
subplot(3,3,7); hold on;
cI = [2 1]; % contrast is reversed
clear pp;
for i = 1:2
    y = doverlap(:,adaptT,cI(i));
    pp(i,:) = squeeze(prms(mi,adaptT,cI(i),:));
    tr(i) = thresh(mi,adaptT,cI(i));
    sl(i) = slope(mi,adaptT,cI(i));
    plot(x,y,'.','color',cols_lite(i,:),'markersize',10);
    plot(xf,mdl(pp(i,:),xf),'color',cols(i,:),'linewidth',1);
    plot([tr(i) tr(i)],[.5 mdl(pp(i,:),tr(i))],...
         '--','color',cols(i,:),'linewidth',1);
end
xlabel('Target Mean');
ylabel('Discriminability');
ylim([.5,1]);xlim([min(muTlow),max(muTlow)])
plotPrefs;

% plot threshold
subplot(3,3,3);
b = bar([1 2],tr,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('Threshold');
plotPrefs; ylim([0 2.5]);

% plot slope
subplot(3,3,6);
b = bar([1 2],sl,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('Slope');
plotPrefs;  ylim([0 .15]);

% plot tau
subplot(3,3,9); hold on;
b = bar([1 2],tau,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('\tau (steps)');
plotPrefs; ylim([0 15]);

saveFigPDF(f12,[350 300],fullfile(plot_dir,'_model_summary.pdf'));



%% summary
mu_lo = 1.5;
mu_hi = mu_lo;
cols = [0 0 1; 1 0 0];
cols_lite = [.8 .8 1; 1 .8 .8];

% extract each time course
dd(:,1) = [rr(muTlow == mu_lo).buffer.Dkl(95:100,1);...
          rr(muTlow == mu_lo).buffer.Dkl(1:50,1)];
dd(:,2) = [rr(muTlow == mu_hi).buffer.Dkl(45:50,1);...
          rr(muTlow == mu_hi).buffer.Dkl(51:100,1)];
kk(:,1) = mean([rr(muTlow == mu_lo).buffer.k(95:100,:);...
                rr(muTlow == mu_lo).buffer.k(1:50,:)],2);
kk(:,2) = mean([rr(muTlow == mu_hi).buffer.k(45:50,:);...
                rr(muTlow == mu_hi).buffer.k(51:100,:)],2);

% time sampling, etc
t = -5:50;
tI = t>=0;
p0 = [.7 .1 4];
yy = dd;

% plot discriminability and fit
f12 = figure(12); clf;
subplot(3,3,[1 2]); hold on
for i = 1:2
    [prm(:,i),emdl,tau(i)] = fitExpGrid(t(tI),movmean(yy(tI,i),2),p0);
    plot(t,yy(:,i),'.','markersize',10,'color',cols_lite(i,:));
    plot(t(tI),emdl(prm(:,i),t(tI)),'color',cols(i,:),'linewidth',1);
end
axis tight; ylim([.5 1]);
plot([0 0],ylim,'k--'); plotPrefs;
ylabel('Discriminability'); xlabel('Time (steps)');

% plot gain
subplot(3,3,[4 5]); hold on
for i = 1:2
    plot(t,kk(:,i),'-','color',cols(i,:));
end
axis tight;
plot([0 0],ylim,'k--');
ylabel('Gain');
xlabel('Time (steps)');
plotPrefs;

% plot psychometric function near full adaptation
clear sl pp tr;
adaptT = 26;
mi = 3;
mdl = fitLogGrid();
subplot(3,3,7); hold on;
cI = [2 1]; % contrast is reversed
clear pp;
for i = 1:2
    y = doverlap(:,adaptT,cI(i));
    pp(i,:) = squeeze(prms(mi,adaptT,cI(i),:));
    tr(i) = thresh(mi,adaptT,cI(i));
    sl(i) = slope(mi,adaptT,cI(i));
    plot(x,y,'.','color',cols_lite(i,:),'markersize',10);
    plot(xf,mdl(pp(i,:),xf),'color',cols(i,:),'linewidth',1);
    plot([tr(i) tr(i)],[.5 mdl(pp(i,:),tr(i))],...
         '--','color',cols(i,:),'linewidth',1);
end
xlabel('Target Mean');
ylabel('Discriminability');
ylim([.5,1]);xlim([min(muTlow),max(muTlow)])
plotPrefs;

% plot threshold
subplot(3,3,3);
b = bar([1 2],tr,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('Threshold');
plotPrefs; ylim([0 2.5]);

% plot slope
subplot(3,3,6);
b = bar([1 2],sl,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('Slope');
plotPrefs;  ylim([0 .15]);

% plot tau
subplot(3,3,9); hold on;
b = bar([1 2],tau,'FaceColor','flat');
b.CData(1,:) = cols(1,:);
b.CData(2,:) = cols(2,:);
ylabel('\tau (steps)');
plotPrefs; ylim([0 15]);

saveFigPDF(f12,[350 300],fullfile(plot_dir,'_model_summary_matchedT.pdf'));


%% plot distributions for each time point in the curve
f13 = figure(13); clf; hold on;
cG = [0,120,87]./255;
%iT = [1,5,9];
times = 0:5:26;
iT = 1:length(times);
bins = -.5:(nLevels-.5);
bb = [bins,nLevels+.5,nLevels+1.5]-1;
for i=1:numel(iT)
    subplot(3,length(iT),i+length(iT));hold on;
    hB = histcounts(rr(muTlow == mu_lo).buffer.yB(100+times(iT(i)),:),bins);
    hB = [0,hB./sum(hB),0,0];
    
    hT = histcounts(rr(muTlow == mu_lo).buffer.yT(100+times(iT(i)),:),bins);
    hT = [0,hT./sum(hT),0,0];
    
    stairs(bb,hB,'color',cols(1,:),'linewidth',1)
    stairs(bb,hT,'color',cG,'linewidth',1)
    ylim([0,.6])
    xlabel('spike count k')
    ylabel('P(k)')
    title(sprintf('Sample %d after switch',times(iT(i))));
    plotPrefs;
    
    subplot(3,length(iT),i+2*length(iT));hold on;
    hB = histcounts(rr(muTlow == mu_hi).buffer.yB(50+times(iT(i)),:),bins);
    hB = [0,hB./sum(hB),0,0];
    
    hT = histcounts(rr(muTlow == mu_hi).buffer.yT(50+times(iT(i)),:),bins);
    hT = [0,hT./sum(hT),0,0];
    
    stairs(bb,hB,'color',cols(2,:),'linewidth',1)
    stairs(bb,hT,'color',cG,'linewidth',1)
    ylim([0,.6])
    xlabel('spike count k')
    ylabel('P(k)')
    plotPrefs;
    
end
clite = [.6 .6 1; 1 .6 .6];
subplot(3,length(iT),[1:3]); hold on;
for i = 1:2
    plot(t,yy(:,i),'-','color',clite(i,:),'linewidth',1);
    plot(times,yy(ismember(t,times),i),'o','markersize',10,...
         'color',cols(i,:));
end
axis tight; ylim([.5 1]);
plot([0 0],ylim,'k--'); plotPrefs;
ylabel('Discriminability'); xlabel('Time (steps)');
saveFigPDF(f13,[900 450],fullfile(plot_dir,'target_noise_dist_dt.pdf'));





%% subplot for each target mean over time inverted and normalized


%  %% ROC computation
%  addpath(genpath('~/chris-lab/code_general/'));
%  for i = 1:length(rr)
%      signalD = rr(i).buffer.yT;
%      noiseD = rr(i).buffer.yB;
%      for t = 1:size(signalD,1)
%          [hits,fas,Dauc(i,t)] = computeROC(noiseD(t,:),signalD(t,:));
%      end
%  end
%  
%  
%  f11 = figure(11); clf; hold on;
%  
%  cmap(:,:,1) = cBL;
%  cmap(:,:,2) = cRL;
%  
%  for i = 1:numel(muTlow)
%      
%      subplot(4,4,i); hold on;
%      clear d;
%      d(1,:) = [Dauc(i,95:100) Dauc(i,1:50)];
%      d(2,:) = [Dauc(i,45:50) Dauc(i,51:100)];
%     
%      plot([0 0],[.5 .9],'k--');
%      title(sprintf('Target mu = %3.2f',muTlow(i)));
%      
%      for j = 1:size(d,1)
%          
%          % plot
%          clear DD;
%          DD = d(j,7:end);
%          
%          if j == 2
%              DD = DD;
%          end
%          %DD = normalize(DD,'range');
%          
%          plot(1:50,DD,'linewidth',2,'color',cmap(i+2,:,j));
%          
%      end
%      plotPrefs;
%          
%  end
%  saveFigPDF(f11,[800 700],'./_plots/_discriminabilityDTperMu_norminv.pdf');



%  for m = 1:length(metrics)
%      
%      fig(m) = figure(m); clf;
%      
%      % plot time markers
%      subplot(5,2,1); hold on;
%      for i = 1:length(tI)
%          plot([times(tI(i)) times(tI(i))], [.5 1],'color',tcolors{1}(i,:),...
%               'linewidth',1);
%          plot([times(tI(i)) times(tI(i))], [1.5 2],'color',tcolors{2}(i,:),...
%               'linewidth',1);
%      end
%      xlabel('Time (samples)'); xlim([-1 26]); plotPrefs;
%      
%      % plot volume markers
%      subplot(5,2,2); hold on;
%      for i = 1:length(muTlow)
%          plot([muTlow(i) muTlow(i)], [.5 1],'color',tcolors{1}(end,:),...
%               'linewidth',mmarks(i));
%          plot([muTlow(i) muTlow(i)], [1.5 2],'color',tcolors{2}(end,:),...
%               'linewidth',mmarks(i));
%      end
%      xlabel('Volume'); xlim([-.25 3.25]); plotPrefs;
%      
%      
%      % extract measure
%      eval(sprintf('D = %s;',metrics{m}));
%      
%      % limits
%      ylims = [floor(min(D(:))*10)/10 ceil(max(D(:))*10)/10];
%      
%      % plot psych curves
%      subplot(5,2,[3 4 5 6]); hold on
%      for i = 1:length(tI)
%          plot(muTlow,squeeze(D(:,tI(i),1)),'color',tcolors{1}(i,:));
%      end
%      for i = 1:length(tI)
%          plot(muTlow,squeeze(D(:,tI(i),2)),'color',tcolors{2}(i,:));
%      end
%      xlabel('Volume'); ylabel(sprintf('D_{%s}',metrics{m})); ylim(ylims);
%      plotPrefs;
%      
%      % plot time curves
%      s = subplot(5,2,[7 9]); hold on
%      for i = 1:length(muTlow)
%          plot(times(tI),squeeze(D(i,tI,1)),'color',tcolors{1}(end,:),...
%               'linewidth',mmarks(i));
%          plot(times(tI),squeeze(D(i,tI,2)),'color',tcolors{2}(end,:),...
%               'linewidth',mmarks(i));
%      end
%      xlabel('Time'); ylabel(sprintf('D_{%s}',metrics{m})); ylim(ylims);
%      %s.Position(3) = s.Position(3) / 2;
%      plotPrefs;
%      
%      subplot(5,2,[8]); hold on;
%      plot(g(:,1),squeeze(sense(m,:,1)),'r');
%      plot(g(:,2),squeeze(sense(m,:,2)),'b');
%      xlabel('Gain'); ylabel('Sensitivity'); plotPrefs;
%      
%      subplot(5,2,[10]); hold on;
%      plot(g(:,1),squeeze(slope(m,:,1)),'r');
%      plot(g(:,2),squeeze(slope(m,:,2)),'b');
%      xlabel('Gain'); ylabel('Slope'); plotPrefs;
%      
%      fn = sprintf('./_plots/model_gain_%s.pdf',metrics{m});
%      saveFigPDF(fig(m),[550 1000],fn);
%      
%  
%  end



        
    
    

    
    





%% plot dependence on target (new figs)

%  %generate colormap
%  cw = [1,1,1];
%  cr = [1,0,0];
%  cb = [0,0,1];
%  cRL = [];cBL = [];
%  for i=1:15
%      cRL = [cRL;cw+(cr-cw)*i./(numel(muTlow)+3)];
%      cBL = [cBL;cw+(cb-cw)*i./(numel(muTlow)+3)];
%  end



%  %plot discriminability over time
%  f1 = figure(1); clf; hold on; 
%  for i=1:numel(muTlow)
%      subplot(1,2,1);hold on;plot(-5:50,[D(95:100,i);D(1:50,  i)],'linewidth',2,'color',cBL(i+2,:));
%      subplot(1,2,2);hold on;plot(-5:50,[D(45:50, i);D(51:100,i)],'linewidth',2,'color',cRL(i+2,:));
%  end
%  subplot(1,2,1);
%  xlabel('time after switch from high to low');
%  ylabel('discriminability');
%  xlim([-5,50]);ylim([.5,.9])
%  plotPrefs;
%  
%  subplot(1,2,2);
%  xlabel('time after switch from high to low');
%  ylabel('discriminability');
%  xlim([-5,50]);ylim([.5,.9])
%  set(gcf,'color','w','Position',[200 200 1200 400])
%  plotPrefs;
%  saveFigPDF(f1,[900 300],'./_plots/_offCurvesdMu.pdf');





    
    








    
    
    
    
    