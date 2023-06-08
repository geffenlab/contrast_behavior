function t = nanpadTable(T,var,grp,grpcnt)

%% function t = nanpadTable(T,var,grp,grpcnt)
%
% T is input table
% var is the table variable to pad
% grp is the group variable to insert nans over

if nargin < 4
    grpcnt = 1;
end

uvar = unique(T.(var));
ugrp = unique(T.(grp));
tvars = T.Properties.VariableNames;
tvars(contains(tvars,{var,grp})) = [];

% get variable sizes
for i = 1:length(tvars)
    varsize{i} = size(T.(tvars{i}));
end

% for each unique variable
t = table;
cnt = 1;
for i = 1:length(uvar)
    
    % check for each group value
    for j = 1:length(ugrp)
        
        if iscell(uvar)
            vI = contains(T.(var),uvar{i});
        else
            vI = T.(var) == uvar(i);
        end
        
        if iscell(ugrp)
            gI = contains(T.(grp),ugrp{j});
        else
            gI = T.(grp) == ugrp(j);
        end

        if any(contains(tvars,'GroupCount'))
            gcI = T.GroupCount < grpcnt;
        else
            gcI = true(size(vI));
        end
        
        I = vI & gI & ~gcI;
        
        if sum(I) == 1
            t(cnt,:) = T(I,:);
            
        elseif sum(I) == 0
            if iscell(uvar)
                t.(var){cnt} = uvar{i};
            else
                t.(var)(cnt) = uvar(i);
            end
            if iscell(ugrp)
                t.(grp){cnt} = ugrp{j};
            else
                t.(grp)(cnt) = ugrp(j);
            end
            
            for k = 1:length(tvars)
                t.(tvars{k})(cnt,:) = nan(1,varsize{k}(2));
            end
            
        end
        cnt = cnt + 1;
    end
end
    
    
    