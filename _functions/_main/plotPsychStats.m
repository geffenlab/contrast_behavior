function [h,b] = plotPsychStats(mp,stat,grp,grpvals,colors,alphas,linestyle,distcolor,ms,labels,connectTheDots)

if size(colors,1) == 1
    colors = repmat(colors,length(grpvals),1);
end

for i = 1:length(grpvals)
    if iscell(grpvals)
        data{i} = mp.(stat)(mp.(grp)==grpvals{i});
    else
        data{i} = mp.(stat)(mp.(grp)==grpvals(i));
    end
end

if ~exist('labels','var') | isempty(labels)
    labels = {[],[]};
end

hold on;
if connectTheDots
    
    % plot the distribution
    h = plotSpread(data,'distributionColors',distcolor);
    
    % get the xvalues
    for i = 1:numel(h{1})
        xv{i} = get(h{1}(i),'XData')';
    end

    % check for missing data
    n = cellfun(@numel,data);
    uM = unique(mp.mouse);
    if any(diff(n) ~= 0)
        
        % make a matrix that is n mice long, filled with nans, and
        % fill it out with the data we have
        pdat = nan(numel(uM),numel(n));
        pxv = pdat;
        
        % for each group
        for i = 1:numel(grpvals)
            
            try
                % find mice with data
                if iscell(grpvals)
                    I = ismember(uM,mp.mouse(mp.(grp)==grpvals{i}));
                else
                    I = ismember(uM,mp.mouse(mp.(grp)==grpvals(i)));
                end
                
                pdat(I,i) = data{i};
                pxv(I,i) = xv{i};
            catch ME
                rethrow(ME);
                keyboard
            end
            
        end
        
    else
        pxv = cell2mat(xv);
        pdat = cell2mat(data);
        
    end
    
    % plot
    plot(pxv',pdat','k','LineWidth',.5);
    
end

% plot remainder
h = plotSpread(data,'distributionColors',distcolor);
b = barWithError(data,labels,colors);
for i = 1:numel(b)
    b(i).FaceAlpha = alphas(i); 
    b(i).LineWidth = 1;
    b(i).EdgeColor = colors(i,:);
    b(i).LineStyle = linestyle{i};
    
    set(h{1}(i),'Marker','o','MarkerFaceColor','w','MarkerSize',ms-1);
end
hold off;

