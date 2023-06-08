function plotFigs_CFA(res)

sigmaNoise = .01;
nLevels = 16;

% cB = [134,152,255]./255;
% cR = [255,148,156]./255;
cB = [0 0 1];
cB_lite = [.6 .6 1];
cR = [1 0 0];
cR_lite = [1 .6 .6];
cG = [0,120,87]./255;

%plot stimulus distributions
ind = 4;
nrows = 4;
ncols = 3;

f1 = figure(1); clf;

s1 = subplot(nrows,ncols,[1:2]);hold on;
plot(1:51,res.buffer.xB(1:51,ind),'-','color',cB_lite);
plot(51:101,res.buffer.xB(51:101,ind),'-','color',cR_lite);
plot(101:150,res.buffer.xB(101:150,ind),'-','color',cB_lite);
plot(1:50,res.buffer.xT(1:50,ind),'.','color',cB,'MarkerSize',10)
plot(51:100,res.buffer.xT(51:100,ind),'.','color',cR,'MarkerSize',10)
plot(101:150,res.buffer.xT(101:150,ind),'.','color',cB,'MarkerSize',10)
xlabel('time')
ylabel('signal')
yticks([-6,-3,0,3,6])
xticks([1,50,100,150])
ylim([-6,6])
plotPrefs;

%---- plot signal distributions ------%
bins = -6:.1:6;

cBH = histcounts(res.xB(res.theta>1.5),bins);
cBH = cBH./sum(cBH);

cBL = histcounts(res.xB(res.theta<1.5),bins);
cBL = cBL./sum(cBL);

cTH = histcounts(res.xT(res.theta>1.5),bins);
cTH = cTH./sum(cTH);

cTL = histcounts(res.xT(res.theta<1.5),bins);
cTL = cTL./sum(cTL);

s2 = subplot(nrows,ncols,3);hold on;
stairs(cBL,bins(1:end-1),'color',cB_lite,'linewidth',1)
stairs(cBH,bins(1:end-1),'color',cR_lite,'linewidth',1)
stairs(cTL, bins(1:end-1),'color',cB,'linewidth',1.5)
stairs(cTH, bins(1:end-1),'color',cR,'linewidth',1.5)
%set(gca,'xscale','log');
%xlabel('log(PDF)')
xlabel('PDF');
ylabel('signal')
yticks([-6,-3,0,3,6])
%xticks(0:.02:.08)
plotPrefs;
%s2.Position(3) = s1.Position(3) / 2;


%----- plot firing rate ---------%

subplot(nrows,ncols,[4:5]);hold on
nT = 3;
plot(1:50,mean(res.buffer.yB(1:50,ind:(ind+nT-1)),2),'color',cB);
plot(101:150,mean(res.buffer.yB(101:150,ind:(ind+nT-1)),2),'color',cB);
plot(51:100,mean(res.buffer.yB(51:100,ind:(ind+nT-1)),2),'color',cR);
plot(mean(res.buffer.yT(:,ind:(ind+nT-1)),2),'color',cG);
xlabel('time')
ylabel('firing rate')
xticks([1,50,100,150])
plotPrefs;


%----- plot discriminability ---------%

subplot(nrows,ncols,[7:8]);hold on
plot(1:150,res.buffer.Dkl(:,1),'--','color',[.6,.6,.6]);
plot(1:50,res.buffer.Dkl(1:50,1),'-o','color',cB);
plot(51:100,res.buffer.Dkl(51:100,1),'-o','color',cR);
plot(101:150,res.buffer.Dkl(101:150,1),'-o','color',cB);
xlabel('time')
ylabel('discriminability')
xticks([1,50,100,150])
plotPrefs;


%----- plot theta estimate ---------%
subplot(nrows,ncols,[10:11]); hold on;
plot(1:150,res.buffer.theta(:,2),'--','color',[.6 .6 .6]);
plot(1:50,res.buffer.theta(1:50,2),'-','color',cB,'linewidth',1);
plot(51:100,res.buffer.theta(51:100,2),'-','color',cR,'linewidth',1);
plot(101:150,res.buffer.theta(101:150,2),'-','color',cB,'linewidth',1);
plot(1:150,mean(res.buffer.thetahat,2),'k','linewidth',1);
xlabel('time');
ylabel('theta');

saveFigPDF(f1,[550 800],'./_plots/_efficient_coding_model.pdf');



%plot nonlinearities
f2 = figure(22);hold on;
times = 1:1:25;
nT = numel(times);
cRmax = [255,208,216]./255;
cBmax = [184,202,255]./255;

cRmin = [160,0,0]./255;
cBmin = [0,0,160]./255;
for i=1:nT
    cmapR(i,1:3) = cRmin+(i-1)*(cRmax-cRmin)./(nT-1);
    cmapB(i,1:3) = cBmin+(i-1)*(cBmax-cBmin)./(nT-1);
end
cmapR = flipud(cmapR);
cmapB = flipud(cmapB);

x0 = -10:.01:10;
for i=1:nT
    %high variance state
    t = times(i);
    
    subplot(2,2,1);hold on;
    [~,rV] = encodeWithDiscrete(x0, res.buffer.k(50+t,ind), res.buffer.x0(50+t,ind), sigmaNoise, nLevels);
    plot(x0,rV,'color',cmapR(i,:),'linewidth',1)
    xlabel('signal')
    ylabel('spike count')
    yticks(linspace(0,1,10))
    yticklabels({'0','1','2','3','4','5','6','7','8','9'})
    xticks(-10:5:10)
    plotPrefs;

    
    % low variance state
    
    subplot(2,2,2);hold on;
    [~,rV] = encodeWithDiscrete(x0, res.buffer.k(100+t,ind), res.buffer.x0(100+t,ind), sigmaNoise, nLevels);
    plot(x0,rV,'color',cmapB(i,:),'linewidth',1)
    yticks(linspace(0,1,10))
    yticklabels({'0','1','2','3','4','5','6','7','8','9'})
    xticks(-10:5:10)
    xlabel('signal')
    ylabel('spike count')
    plotPrefs;

end


bins = -10:.1:10;

cBH = histcounts(res.xB(res.theta>1.5),bins);
cBH = cBH./sum(cBH);

cBL = histcounts(res.xB(res.theta<1.5),bins);
cBL = cBL./sum(cBL);

subplot(2,2,3);hold on;
stairs(bins(1:end-1),cBL,'color',[.6,.6,.6])
stairs(bins(1:end-1),cBH,'color',cR,'linewidth',1)
ylabel('PDF')
xlabel('signal')
xticks([-10,-5,0,5,10])
yticks(0:.02:.04)
plotPrefs;



subplot(2,2,4);hold on;
stairs(bins(1:end-1),cBH,'color',[.6,.6,.6])
stairs(bins(1:end-1),cBL,'color',cB,'linewidth',1)
ylabel('PDF')
xlabel('signal')
xticks([-10,-5,0,5,10])
yticks(0:.02:.04)
plotPrefs;

saveFigPDF(f2,[400 300],'./_plots/nl_adapt.pdf')



%plot histograms

f3 = figure(33); clf; hold on;
%iT = [1,5,9];
times = 1:3:25;
iT = 1:length(times);
bins = -.5:(nLevels-.5);
bb = [bins,nLevels+.5,nLevels+1.5]-1;
for i=1:numel(iT)
    subplot(2,length(iT),i);hold on;
    hB = histcounts(res.buffer.yB(50+times(iT(i)),:),bins);
    hB = [0,hB./sum(hB),0,0];
    
    hT = histcounts(res.buffer.yT(50+times(iT(i)),:),bins);
    hT = [0,hT./sum(hT),0,0];
    
    stairs(bb,hB,'color',cR,'linewidth',1)
    stairs(bb,hT,'color',cG,'linewidth',1)
    ylim([0,.6])
    xlabel('spike count k')
    ylabel('P(k)')
    title(sprintf('Sample %d after switch',times(iT(i))));
    plotPrefs;
    
    subplot(2,length(iT),i+length(iT));hold on;
    hB = histcounts(res.buffer.yB(100+times(iT(i)),:),bins);
    hB = [0,hB./sum(hB),0,0];
    
    hT = histcounts(res.buffer.yT(100+times(iT(i)),:),bins);
    hT = [0,hT./sum(hT),0,0];
    
    stairs(bb,hB,'color',cB,'linewidth',1)
    stairs(bb,hT,'color',cG,'linewidth',1)
    ylim([0,.6])
    xlabel('spike count k')
    ylabel('P(k)')
    plotPrefs;
    
end
saveFigPDF(f3,[1600 300],'./_plots/target_noise_dist_dt.pdf');


%plot discriminability
f4 = figure(44); clf; hold on;

subplot(2,1,1);hold on;
plot((res.buffer.Dkl(1:50,1)),'-o','color',cB);
plot((res.buffer.Dkl(51:100,1)),'-o','color',cR);
for i=1:numel(iT)
    plot(times(iT(i)),res.buffer.Dkl(50+times(iT(i)),1),'o','markersize',16,'linewidth',2,'color',cmapR(iT(i),:))
    plot(times(iT(i)),res.buffer.Dkl(100+times(iT(i)),1),'o','markersize',16,'linewidth',2,'color',cmapB(iT(i),:))
end
xlabel('time after switch')
ylabel('discriminability')
plotPrefs;

subplot(2,1,2);hold on;
for i=1:numel(times)
    plot([times(i),times(i)],[1,2],'-','linewidth',2,'color',cmapB(i,:))
    plot([times(i),times(i)],[3,4],'-','linewidth',2,'color',cmapR(i,:))
end
xlim([0,50])
ylim([0,5])
xlabel('time after switch')
plotPrefs;

saveFigPDF(f4,[300 300],'./_plots/discriminability_dt.pdf');

