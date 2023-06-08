function plotPvalSmart(p,x,y,value)

hold on;
units = 'normalized';
if numel(x) > 1
    % for line plots
    xl = x;
    x = mean(x);
    YL = ylim;
    yl = [y y] + range(YL)*.1;
    y = y + range(YL)*.15;
    if p < .05
        plot(xl,yl,'k','LineWidth',1);
    end
    units = 'data';
end

[psym,pval] = pvalStr(p);
if exist('value','var') & ~isempty(value)
    if strcmp(value,'pval')
        str = sprintf('p=%s',pval);
    elseif strcmp(value,'psym')
        str = sprintf('%s',psym);
    end
else
    str = sprintf('p=%s %s',pval,psym);
end

h = text(x,y,str,'units',units,'horizontalAlignment','center');
hold off;

