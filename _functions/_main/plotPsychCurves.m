function plotPsychCurves(t,I,color,linestyle,markerfill,ms,thresh2use,indcurves)

% fit
x = t.mean_vols(I,:); xf = linspace(min(x(:)),max(x(:)),100);
y = t.mean_pc(I,:);
[prms,mdl,thresh,~,~,~,thresh75] = fitLogGrid(x(:),y(:),[],[],[],.75);

if strcmp(thresh2use,'thresh75')
    thresh = thresh75;
end

% plot individual curves
if ~exist('indcurves','var') | indcurves
    plot(x',y','Color',color + ((color == 0) * .7),'LineStyle', ...
         linestyle);
else
    errorBars(mean(x,1,'omitnan'),y,color + ((color == 0) * .7),'sem',[],[],[],...
              'o','LineWidth',1,'MarkerSize',ms, ...
              'MarkerFaceColor',markerfill);
end

%  % plot data
%  errorBars(mean(x,1,'omitnan'),y,color + ((color == 0) * .7),'sem',[],[],[],...
%            'o','LineWidth',1,'MarkerSize',ms, ...
%            'MarkerFaceColor',markerfill);

% fits
plot(xf,mdl(prms,xf),'Color',color,'LineWidth',1,'LineStyle',linestyle);
plot([thresh thresh],[.5 mdl(prms,thresh)],'--','Color',color,'LineStyle',linestyle)
%plot(thresh,mdl(prms,thresh),'k.','MarkerSize',ms*2);