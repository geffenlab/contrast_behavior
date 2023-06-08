function [t,uOff] = buildOffsetBehaviorRes(dat,fileName)

t = table;
uOff = [];

if ~exist(fileName,'file')
    % extract offset sessions
    t = dat(contains(dat.sesstype,'offset'),:);
    t = cleanTable(t,{'vols','volP'});


    %% reformat for all unique offset times
    % variables to update
    var2fix = {'offs','rate','dp','pc','nresp','ntrials'};

    % unique offsets
    uOff = unique(round([t.offs{:}]-3,3));

    % for each session, remap offsets to the unique offsets
    for i = 1:size(t,1)
        
        fprintf('Sess %d/%d... ',i,size(t,1)); tic;
        
        % get current offset as an index into uOffs
        offI = ismember(uOff,round(t.offs{i}-3,3));
        
        % for each variable to fix, make nans then add indexes
        for j = 1:length(var2fix)
            
            val = t.(var2fix{j}){i};
            
            % for each row in the variable, duplicate it to new vars
            for k = 1:size(val,1)
                
                % get values
                newvar = sprintf('%s_%d',var2fix{j},(k-1));
                rval = t.(var2fix{j}){i}(k,:);
                
                % assign remapped measures to new variables
                t.(newvar)(i,:) = nan(1,length(uOff));
                t.(newvar)(i,offI) = rval;
                
            end
        end
        
        % fit tau for each session
        x = t.offs{i};
        y = t.pc{i}(1,:);
        [params,mdl,TAU] = fitExpGrid(x,y);
        t.tau(i) = TAU;
        t.exp_fit(i,:) = params';
        toc;
        
    end
    t.offs_0 = t.offs_0 - 3;
    
    save(fileName,'t','uOff');

else
    
    load(fileName);
    
end