function plotMultiComp(stats,y)

% plots pvalues and lines for multiple comparison tests, it assumes
% that each group is plotted at 1,2,...,n for n groups
t = multcompare(stats,'display','off');

% draw significant stats
pcorr = .05 / size(t,1);
I = find(t(:,end)<pcorr);
YL = ylim;
hold on;
for i = 1:length(I)
    plot(t(I(i),1:2),[y y] + (i)*range(YL)*.1,'k','LineWidth',1);
    psym = pvalStr(t(I(i),end));
    text(mean(t(I(i),1:2)), y + (i)*range(YL)*.125,...
         sprintf(psym),'horizontalAlignment','center');
end
hold off;