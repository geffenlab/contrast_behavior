function [m,e,dat] = combindstats(data,ind);

%% function [m,e,dat] = combindstats(data,ind);
%
% this function combines columns of data according to the values in
% ind. ind should be a cell array, with each cell containing a
% vector indicating the columns to be averaged for that index
%
% returns:
%  m - mean over not-nan rows for each index
%  s - sem over not-nan rows for each index
%  dat - data grouped by ind
%
% example: (if data is a 10 column array, the following index will
% average each pair of columns):
% 
% data = rand(100,10);
% clear ind;
% ind{1} = [1 2];
% ind{2} = [3 4];
% ind{3} = [5 6];
% ind{4} = [7 8];
% ind{5} = [9 10];
% [m,e,dat] = comdindstats(data,ind);


for i = 1:length(ind)
    colI = ind{i};
    tmp = data(:,colI);
    m(i) = mean(tmp(:),'omitnan');
    e(i) = std(tmp(:),'omitnan') ./ sqrt(sum(~isnan(tmp(:))));
    dat{i} = tmp;
end