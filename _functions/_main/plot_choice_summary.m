function stat = plot_choice_summary(res_choice,res_ln,sessionData,ops)


% align results for choice and LN recordings
[C,ia,ib] = intersect({res_choice.cellID},{res_ln.cellID});
rc = res_choice(ia);
rl = res_ln(ib);

% gain
gain = cat(3,rl.NRC_ahat);
gain = squeeze(gain(:,3,:))';

% choice probability
cp = cat(3,rc.cp);
mcp = squeeze(mean(cp,1,'omitnan'))';
cp_sig = cat(3,rc.cp_sig);
ncp_sig = squeeze(sum(cp_sig,2))';

% contrast
contrast = contains({rc.contrast},'lohi')';

% noise ratio
NR = squeeze(median(cat(3,rl.noiseRatio),1,'omitnan'))';
include = all(NR > 0 & NR < 100,2) & sum(ncp_sig,2) > 20 & ...
          contrast == 0;

% gain and cp correlations
[targ_corr,targ_pv] = corr(gain(include,2),mcp(include,1:end-1));
[adapt_corr,adapt_pv] = corr(gain(include,1),mcp(include,1:end-1));

% plot
figure(123); clf; hold on;
ts = ops.timeCent(1:end-1);
bf_alpha = 0.05 / 20;
plot(ts(targ_pv<bf_alpha),targ_corr(targ_pv<bf_alpha),...
     'r-')
plot(ts,targ_corr,'k');
plot(ts(adapt_pv<bf_alpha),adapt_corr(adapt_pv<bf_alpha),...
     'r-');
plot(ts,adapt_corr,'k--');




% reaction time analysis
sI = contains({sessionData.session.task},'psychometric');
s = sessionData.session(sI);
b = sessionData.behavior(sI);

for i = 1:numel(b)
    
    % for each session, get reaction times for each volume
    rtI = b(i).RT > .05 & b(i).RT < 1;
    rt = b(i).RT(rtI);
    tt = b(i).lickTT(rtI);
    snrI = [0 1 2 3 4 5 6];
    currentSNRS = ismember(snrI,tt);
    RT(i,:) = nan(size(currentSNRS));
    snr(i,:) = [-inf s(i).stimInfo.targetDBShift];
    RT(i,currentSNRS) = grpstats(rt,tt);
    snri(i,:) = snrI;
    contrast(i) = contains(s(i).cond,'lohi');
    
end
cc = repmat(contrast',1,7);
grpRT_mean = grpstats(RT(:),{snr(:),cc(:)},'mean');
grpRT_sem = grpstats(RT(:),{snr(:),cc(:)},'sem');
grps = grpstats([snr(:),cc(:)],{snr(:),cc(:)});
grpSNR = grpstats(snr(:),{snr(:),cc(:)});      

figure;
% group RTs by contrast and volume
color = {'b','r'};
hold on
for i = 1:2
    x = grpSNR(grps(:,2)==(i-1));
    ym = grpRT_mean(grps(:,2)==(i-1));
    ye = grpRT_sem(grps(:,2)==(i-1));
    errorbar(-10,ym(1),ye(1),...
             color{i},'linewidth',1,'marker','o');
    errorbar(x(2:end),ym(2:end),ye(2:end),...
             color{i},'linewidth',1,'marker','o');
    
end
xlim([-12 27]); 
ylabel('RT (s)'); xlabel('Target Volume (mean dB SNR)');
plotPrefs;