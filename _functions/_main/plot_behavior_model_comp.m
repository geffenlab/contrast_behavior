function fh = plot_behavior_model_comp(behavior,model,fun2use)

if ~exist('fun2use','var') | isempty(fun2use)
    fun2use = 1;
end

% percent change function
if fun2use == 1
    pcfun = @(x) ((x(2) - x(1)) ./ x(1));
    lim = [-1 1.25];
    ylab = '% change';
else
    pcfun = @(x) ((x(2) - x(1)) ./ sum(x(:)));
    lim = [-.55 .55];
    ylab = 'CMI';
end

% bootstrap errorbars on the change in behavior
its = 2000;

% behavioral parameters
thr = cat(2,behavior(5).xx{:});
slp = cat(2,behavior(6).xx{:});
tau = behavior(end).xx;

for i = 1:its
    
    I = randsample(1:size(thr,1),size(thr,1),true);
    
    % thresholds
    tmp_thr = mean(thr(I,:),'omitnan');
    thrd(i) = pcfun(tmp_thr);
    
    % slopes
    tmp_slp = mean(slp(I,:),'omitnan');
    slpd(i) = pcfun(tmp_slp);
    
    I = randsample(1:size(tau,1),size(tau,1),true);
    
    % taus
    tmp_tau = mean(tau(I,:),'omitnan');
    taud(i) = pcfun(tmp_tau);

end

fh = figure(898); clf
subplot(1,3,1); hold on;
bar(2,pcfun(model.threshold),'facecolor',[.5 .5 .5])
bar(1,median(thrd,'omitnan'),'facecolor','k');
plot([1 1],prctile(thrd,[5 95]),'k')
ylabel(sprintf('Threshold (%s)',ylab)); title('Threshold');
set(gca,'xtick',[1 2],'xticklabels',{'Data','Model'})
ylim(lim);
plotPrefs;

subplot(1,3,2); hold on;
bar(2,pcfun(model.slope),'facecolor',[.5 .5 .5])
bar(1,median(slpd,'omitnan'),'facecolor','k');
plot([1 1],prctile(slpd,[5 95]),'k')
ylabel(sprintf('Slope (%s)',ylab)); title('Slope');
set(gca,'xtick',[1 2],'xticklabels',{'Data','Model'})
ylim(lim);
plotPrefs;


subplot(1,3,3); hold on;
bar(2,pcfun(model.tau),'facecolor',[.5 .5 .5])
bar(1,median(taud,'omitnan'),'facecolor','k');
plot([1 1],prctile(taud,[5 95]),'k')
ylabel(sprintf('tau (%s)',ylab)); title('Adaptation Time');
set(gca,'xtick',[1 2],'xticklabels',{'Data','Model'})
ylim(lim);
plotPrefs;
