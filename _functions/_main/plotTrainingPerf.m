function plotTrainingPerf(t,mouseList,colors)

% for each mouse, find the date of the first high contrast and low
% contrast session, then make a vector of performance per session
for i = 1:length(mouseList)
    
    % first date and date vector
    mI = contains(t.mouse,mouseList{i});
    if size(mI,1) == 1
        mI = mI';
    end
    
    % get mean contrast for the 20 sessions
    cc = t.contrast(mI);
    firstHighContrast(i) = round(mean(cc(1:20)));
    
    for j = 1:2
        
        I = t.contrast == j-1 & mI;
        firstDate{i,j} = t.date(find(I,1,'first'));
        expDays{i,j} = round(days(t.date(I) - firstDate{i,j}));
        perfDays{i,j} = t.mean_pc(I);
        n(i,j) = sum(I);
        
    end
    
    
end

% make this padded with nans to have better plots
dayLabel = [0 max(n(:))];
perfMat = nan(max(n(:))+1,length(mouseList),2);
for i = 1:length(mouseList)
    for j = 1:2
        perfMat([expDays{i,j}+1]',i,j) = perfDays{i,j};
    end
end
perfMat(perfMat==0) = nan;


%  firstHighContrast = days(vertcat(firstDate{:,1}) - vertcat(firstDate{:,2})) ...
%      > 0;
%  
%  % mouse CA049 had a few early high contrast days, manually
%  % switch them to low contrast starts
%  firstHighContrast(contains(mouseList,'CA049')) = ...
%      ~firstHighContrast(contains(mouseList,'CA049'));


% plot individual lines for each mouse conditioned on whether that
% mouse started low or high contrast first
hold on
for i = 1:2
    I = find(firstHighContrast==(i-1));
    
    for j = 1:length(I)
        if i == 2 & contains(mouseList{I(j)},'CA121')
            plot(expDays{I(j),i},movmean(perfDays{I(j),i},7,'omitnan'),...
                 '-.','Color',[colors(i,:) .25]);
        else
            plot(expDays{I(j),i},movmean(perfDays{I(j),i},7,'omitnan'),...
                 '-','Color',[colors(i,:) .25]);
        end
    end
    
    ph(i) = plot(movmean(nanmean(perfMat(:,I,i),2),7,'omitnan'),...
                 'Color',colors(i,:),'LineWidth',1);
    %phh(i) = plot(movmean(nanmean(perfMat(:,I,i),2),7,'omitnan'),...
%             '.','Color',colors(i,:),'LineWidth',1);
    text(.8,.1+(.1*(i-1)),sprintf('n=%d',length(I)),...
         'units','normalized','color',colors(i,:));
end
ch = plot(xlim,[.5 .5],'k'); axis tight;
uistack(ph,'top');
%uistack(phh,'top');
xlabel('Time (days from task exposure)');
ylabel('Percent Correct'); plotPrefs;
hold off;
        
        

if 1 == 2
    % plot each session as a dot, with median overlaid
    hold on;
    for i = 1:2
        I = firstHighContrast==(i-1);
        plot(perfMat(:,I,i),'.','Color',colors(i,:)+(colors(i,:)==0)*.8);
        ph(i) = plot(movmean(nanmedian(perfMat(:,I,i),2),7,'omitnan'),...
                     'Color',colors(i,:),'LineWidth',1);
        text(.8,.1+(.1*(i-1)),sprintf('n=%d',sum(I)),...
             'units','normalized','color',colors(i,:));
    end
    ch = plot(xlim,[.5 .5],'k'); axis tight;
    uistack(ph,'top');
    xlabel('Time (days rel. task exposure)');
    ylabel('Percent Correct'); plotPrefs;
    hold off;
end