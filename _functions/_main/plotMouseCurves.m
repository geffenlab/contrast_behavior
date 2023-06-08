function [mouseSNR,mousePerf,mouseContrast,mouseID] = plotMouseCurves(t,mouseList,mdl,colors,ms)

mouseSNR = []; mousePerf = []; mouseContrast = []; mouseID = [];
for i = 1:numel(mouseList)
    
    subplot(5,5,i); hold on;
    plot([-10 30],[.5 .5],'k');
    ph = [nan nan];
    for j = 1:2
        
        % extract data for this mouse and contrast
        mI = contains(t.mouse,mouseList{i});
        if size(mI,1) == 1
            mI = mI';
        end
        I = t.contrast == j-1 & mI;
        
        % exclude sessions with low performance at high volume
        I = I & t.rate(:,end) > .75;
        
        if sum(I) > 0
            x = round(t(I,:).vols,1);
            xf = linspace(min(x(:)),max(x(:)),100);
            y = t(I,:).pc;
                        
            % save out mean performance per SNR for later
            usnrs = grpstats(x(:),x(:),'mean');
            mouseSNR = [mouseSNR; usnrs];
            mousePerf = [mousePerf; grpstats(y(:),x(:),'mean')];
            mouseContrast = [mouseContrast; repmat(j-1,length(usnrs),1)];
            mouseID = [mouseID; repmat(i,length(usnrs),1)];
            
            % plot it
            plot(x',y','-','color',colors(j,:)+(colors(j,:)==0)*.7,...
                 'MarkerSize',ms)
            
            % plot median of individual session fits
            % mp = median(t.prms(I,:)); mt = mp(1)/mp(2);
            
            % fit to all sessions
            [mp,mdl,mt] = fitLogGrid(x(:),y(:));
            
            % plot fit
            ph(j) = plot(xf,mdl(mp,xf),'color',colors(j,:),'LineWidth',1);
            plot([mt mt],[.5 mdl(mp,mt)],'color',colors(j,:));
            plotPrefs;
            drawnow;
        end
        
    end
    for j = 1:length(ph)
        if ~isnan(ph(j))
            uistack(ph(j),'top');
        end
    end
    set(gca,'xtick',[-5 0 5 10 15 20 25]);
    title(mouseList{i});
    axis tight; ylim([.4 1]);
    
    if i == numel(mouseList)
        xlabel('Target Volume (dB SNR)');
        ylabel('Percent Correct');
    end
end