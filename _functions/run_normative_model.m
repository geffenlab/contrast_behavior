function r = run_normative_model

%% initialize 
fprintf('--------------\nRunning normative model simulations...\n');

%define mean of target distribution
muTlow = 0:.25:3;

% contrast widths
sigmaLow = 1;
sigmaHigh = 2;

% target scaling factor
f = .25;

%index of distribution used in previous plots
indTarget = find(muTlow == 1.5);

% visualize pdfs
x = linspace(min(muTlow)-2*sigmaHigh,max(muTlow)+2*sigmaHigh,1000);
for i = 1:length(muTlow)
    subplot(2,1,1); hold on;
    plot(x,normpdf(x,muTlow(i),sigmaLow*f),'b');
    axis tight;
    subplot(2,1,2); hold on;
    plot(x,normpdf(x,muTlow(i),sigmaHigh*f),'r');
    axis tight;
end



%% generate results for all values of muT
resFile = './_data/_res_normative_model.mat';
%resFile = 'discriminability.mat';


%-----generate results for muTlow = 1.5------%
[res,p] = varianceSims(muTlow(indTarget),sigmaLow,sigmaHigh,f);

%% plot figures used in Cosyne abstract
plotFigs_CFA(res); drawnow;



%% generate res for all target means
if ~exist(resFile)
    fprintf('Simulating different target volumes...\n');
    for i = 1:length(muTlow)
        rr(i) = varianceSims(muTlow(i),sigmaLow,sigmaHigh,f);
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
    ylim([.6,.9]);xlim([min(muTlow),max(muTlow)])
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
saveFigPDF(f2,[900 400],'./_plots/_psychCurvesDT.pdf');


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
        
    end
    plotPrefs;
        
end
saveFigPDF(f11,[800 700],'./_plots/_discriminabilityDTperMu.pdf');






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
saveFigPDF(f12,[450 275],'./_plots/deltaGain.pdf');



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
axis tight; ylim([.65 .85]);
plot([0 0],ylim,'k--'); plotPrefs;
ylabel('Discriminability'); xlabel('Time (steps)');

% plot gain
subplot(3,3,[4 5]); hold on
yyaxis left;
plot(t,kk(:,1),'-','color',cols(1,:));
axis tight;
ylabel('Gain_{L}');
yyaxis right
plot(t,kk(:,2),'-','color',cols(2,:));
set(gca,'ydir','reverse')
ylabel('Gain_{H} (inverted)');
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
ylim([.6,.9]);xlim([min(muTlow),max(muTlow)])
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

% save res to output
r.threshold = tr;
r.slope = sl;
r.tau = tau;
r.level_index = mi;
r.time_index = adaptT;


saveFigPDF(f12,[350 300],'./_plots/_model_summary.pdf');


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
axis tight; ylim([.65 .85]);
plot([0 0],ylim,'k--'); plotPrefs;
ylabel('Discriminability'); xlabel('Time (steps)');
saveFigPDF(f13,[900 450],'./_plots/target_noise_dist_dt.pdf');




    
    








    
    
    
    
    