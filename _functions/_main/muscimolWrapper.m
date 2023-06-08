function [res,ops] = muscimolWrapper(d,ops)

ncells = size(d.cellInfo,1);

% process the events a little bit
ev = d.block(1).stimOn;
ev = round(ev - ev(1),3);
stimI = round(mod(ev,5)) == 0;
targI = round(mod(ev,5)) ~= 0;

% stim info
stimInfo = d.block(1).stimInfo.p(1);
stimInfo.order = d.block(1).stimInfo.order;

% labels
blockStr{1} = ['pre' ops.condShort];
blockStr{2} = ['post' ops.condShort];
volLabels = num2str([-inf stimInfo.targetDBShift]');
volTicks = [-5 stimInfo.targetDBShift];

% plot settings and colors
cc = {'b','r'};
rgb = [0 0 1; 1 0 0];
cmap{1} = [.5 .5 .5; linspace(.75,0,6)' linspace(.75,0,6)' ones(6,1)];
cmap{2} = [.5 .5 .5; ones(6,1) linspace(.75,0,6)' linspace(.75,0,6)'];
ms = 4;

for c = 1:ncells
        
    if ~isempty(ops.cellPlot)
        % subplots
        f1 = figure(1); clf;
         
        % make axes
        for b = 1:2
            ax(1+5*(b-1)) = subplot(5,4,[1 2 5 6]+2*(b-1));
            ax(2+5*(b-1)) = subplot(5,4,[9 10]+2*(b-1));
            ax(3+5*(b-1)) = subplot(5,4,[13 14]+2*(b-1));
            ax(4+5*(b-1)) = subplot(5,4,17+2*(b-1));
            ax(5+5*(b-1)) = subplot(5,4,18+2*(b-1));
        end
    end
        
    % for pre and post topical application blocks (block 1 is pre)
    clear YL;
    for b = 1:2
        
        %% compute
        % axis iterator
        axI = 5*(b-1);
        
        % extract spikes
        clustID = d.block(b).clustID(c);
        clustLabel = d.block(b).clustLabel{c};
        spikes = d.block(b).spikes(d.block(b).clust == clustID);
        
        % title
        if b == 1
            tStr = sprintf('%s: cell%d (%s)\npre %s',...
                           ops.mouse,clustID,clustLabel,ops.condLong);
        else
            tStr = sprintf('post %s',ops.condLong);
        end
        
        % make a psth for target events
        tEv = d.block(b).stimOn(stimI) + stimInfo.noiseD;
        [PSTH,raster,trials,PSTHsm] = makePSTH(spikes,tEv, ...
                                               ops.psthWindow,1);
        
        %spkcnt = round(PSTH .* ops.bin);
        %binI = ops.psthWindow >= ops.targetWindow(1) & ...
%       ops.psthWindow <= ops.targetWindow(2);
%spksPerTrial = sum(spkcnt(:,binI),2) ./ diff(ops.targetWindow);
        
        % sort
        [sort,sortI] = sortrows(stimInfo.order);
        [~,spikeI] = ismember(trials,sortI);
        
        % compute averages for each condition
        tI = ops.psthTime > ops.targetWindow(1) & ...
             ops.psthTime < ops.targetWindow(2);
        for i = 1:2
            for j = 1:7
                % mean trace
                condTrace(:,i,j) = mean(PSTHsm(stimInfo.order(:,1)==i& ...
                                               stimInfo.order(:,2)==j-1,:));
                
                % mean over target time index
                condTarget(:,i,j) = mean(squeeze(PSTHsm(stimInfo.order(:,1)==i& ...
                                                        stimInfo.order(:,2)==j-1,tI)),2);
                
                % ROC analysis
                if j == 1
                    noiseD = squeeze(condTarget(:,i,j));
                else
                    signalD = squeeze(condTarget(:,i,j));
                    [auc(i,j-1),~,~,~,aucPct(i,j-1,:)] = ...
                        bootROC(noiseD,signalD);
                    aucSig(i,j-1) = ~(.5 >= aucPct(i,j-1,1) & .5 <= aucPct(i,j-1,2));
                                    
                end
            end
        end
                
        % save res
        res(c).PSTH{b} = PSTH;
        res(c).order = stimInfo.order;
        res(c).targetIndex = tI;
        res(c).condTrace{b} = condTrace;
        res(c).condTarget{b} = condTarget;
        res(c).auc{b} = auc;
        res(c).aucPct{b} = aucPct;
        res(c).aucSig{b} = aucSig;
        
        
        
        
        
        
        
        
        
        
        
        
        
        if ~isempty(ops.cellPlot)
            
            %% plot
            % for each contrast and level
            for i = 1:2
                for j = 1:7
                    axes(ax(1+axI)); hold on;
                    % plot color patch
                    I = [find(all(sort(:,1:2) == [i j-1],2),1,'first')-1 ...
                         find(all(sort(:,1:2) == [i j-1],2),1,'last')];
                    ph = patch(repmat(ops.psthWindow(1),1,4) ./ [1 1 2 2],...
                               [I fliplr(I)],cmap{i}(j,:),'EdgeAlpha', ...
                               0);
                    % and line
                    plot([ops.psthWindow(1) ops.psthWindow(end)],[I(1) I(1)],...
                         'Color',[.75 .75 .75],'LineWidth',.25);
                    if i==2 & j==1
                        plot([ops.psthWindow(1) ops.psthWindow(end)],[I(1) I(1)],...
                             'k','LineWidth',.5);
                    end
                    xlim([ops.psthWindow(1) ops.psthWindow(end)]);
                    
                    % plot PSTH plots
                    if i == 1
                        axes(ax(2+axI)); hold on;
                    else
                        axes(ax(3+axI)); hold on;
                    end
                    plot(ops.psthTime,condTrace(:,i,j),'color', ...
                         cmap{i}(j,:));
                    ylabel('spks/s');
                    axis tight; plotPrefs;
                    xlim([ops.psthWindow(1) ops.psthWindow(end)]);

                end
                
                % mean FR per SNR
                if i == 1
                    axes(ax(4+axI)); hold on;
                else
                    axes(ax(5+axI)); hold on;
                end
                
                plot(volTicks(2:end),squeeze(nanmean(condTarget(:,i,2:end))),...
                     cc{i});
                for j = 1:length(volTicks)
                    eh = errorBars(volTicks(j),squeeze(condTarget(:,i,j)),...
                                   cmap{i}(j,:));
                    eh.Marker = 'o'; eh.MarkerSize = ms;
                    eh.MarkerFaceColor = cmap{i}(j,:);
                end
                set(gca,'xtick',volTicks); set(gca,'xticklabels',volLabels);
                xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
                axis tight;
                YL(b,i,:) = get(gca,'YLim'); plotPrefs;
                xlim([volTicks(1) - mean(diff(volTicks)) ...
                      volTicks(end) + mean(diff(volTicks))]);
                

            end
            
            % plot cleanup
            axes(ax(1+axI)); hold on;
            plot([0 0],[0 length(stimInfo.order)],'k','LineWidth',.5);
            scatter(raster,spikeI,10,'k.');
            title(tStr);
            axis tight; plotPrefs;
        end
        
    end
       
    if ~isempty(ops.cellPlot)
        % match trace plot lims
        tr = [res(c).condTrace{:}];
        ymt = [min(tr(:)) max(tr(:))];
        if ymt(2) == 0; ymt(2) = 1; end;
        set(ax(2),'Ylim',ymt,'XLim',[ops.psthWindow(1) ops.psthWindow(end)]);
        set(ax(3),'Ylim',ymt,'XLim',[ops.psthWindow(1) ops.psthWindow(end)]);
        set(ax(7),'Ylim',ymt,'XLim',[ops.psthWindow(1) ops.psthWindow(end)]);
        set(ax(8),'Ylim',ymt,'XLim',[ops.psthWindow(1) ops.psthWindow(end)]);
        
        % match volume function plots limits
        ym = [min(YL(:)) max(YL(:))];
        if ym(2) == 0; ym(2) = 1; end;
        set(ax(4),'Ylim',ym);
        set(ax(5),'Ylim',ym);
        set(ax(9),'Ylim',ym);
        set(ax(10),'Ylim',ym);
        
        % add zero lines at appropriate ylim
        for i = 2:3
            axes(ax(i)); hold on; plot([0 0],ylim,'k');
            axes(ax(i+5)); hold on; plot([0 0],ylim,'k');
        end
        
        fn = fullfile(ops.cellPlot,...
                      sprintf('%s-cell%d-%s_%s.pdf',...
                              ops.mouse,clustID,clustLabel, ...
                              ops.condShort));
        saveFigPDF(f1,[400 600],fn);
    end
    
end

        
    