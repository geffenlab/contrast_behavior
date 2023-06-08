function buildBehaviorRes(mouseList,filePath,resFile,pcFlag,adjust);

%% builds a data table for each behavioral session

% clear vars
clear fileList fileInd p;

% initialize table
dat = table;

% session types
types = {'noresponse','training','psychometric','offset','opto-psychometric','threshold-opto','switch'};

% for each mouse
for m = 1:length(mouseList)
    
    % mouse directory
    mDir = fullfile(filePath,mouseList{m});
    
    % index data files in that location
    [fileList fileInd prms] = indexDataFiles(mDir);
    fprintf('MOUSE %s\n-----------\n',mouseList{m});
    

    clear s;
    
    % for each session, do analysis
    for i = 1:length(fileList)
        
        fprintf('\t%d/%d - %s ',i,length(fileList), ...
                num2str(fileInd(i,end))); tic;
        
        
        
        %% load session data
        abort = [];
        load(fileList{i});
        
        if isempty(abort)
            abort = zeros(size(resp));
        end
        p = params;
        
        
        
        %% session info
        s(i).mouse = p.IDstr;
        s(i).date = datetime(num2str(fileInd(i,end)),...
                             'InputFormat','yyMMddHHmm');
        s(i).contrast = p.sd(2) > p.sd(1);
        s(i).vols = p.targetDBShift;
        s(i).offs = p.noiseD;
        s(i).noiseD = p.baseNoiseD;
        s(i).mu = p.mu;
        s(i).cwidth = p.sd;
        s(i).stimfile = p.stim;
        s(i).datafile = fileList{i};
        s(i).timeout = p.timeoutD;
        s(i).sesscode = fileInd(i,2);
        s(i).sesstype = types{fileInd(i,2)+1};
        
        % fit the psychometric data
        if length(s(i).vols) == 1 & s(i).sesscode == 2
            s(i).vols = p.attenuation;
        else
            s(i).vols = s(i).vols;
        end
        
        % offset correction (offsets before July 17, 2018 came 25ms early)
        if s(i).date < datetime(2018,07,17)
            % offsets come 25ms early
            s(i).offs = s(i).offs - .025;
            
        end
        
        % reward contingencies
        if ~isfield(p,'rewCont')
            rc = [0 ones(size(s(i).vols))];
        else
            rc = p.rewCont;
        end
        s(i).rewcont = rc;
        
        % volume presentation probability
        if ~isfield(p,'dbP')
            dbp = [.5 .5];
        else
            dbp = p.dbP;
        end
        s(i).volP = dbp;
        
        % offset probability
        if ~isfield(p,'offsetP')
            offp = [.25 .25 .25 .25];
        else
            offp = p.offsetP;
        end
        s(i).offP = offp;
        
        % get file notes
        clear note label;
        if isfield(p,'note') & ~isempty(p.note)
            s(i).note = p.note;
        else
            s(i).note = 'none';
        end
        
        % stimulus label
        if isfield(p,'stimLabel')
            if iscell(p.stimLabel)
                s(i).label = 'none';
            else
                s(i).label = p.stimLabel;
            end
        else
            s(i).label = 'none';
        end

        
        
        

        % compute performance for each trial type (ignore the last
        % column of tt for now -- this is the noise background
        % pattern, and may be incorrect for some mice)
        [rate, dp, pc, nresp, ntrials] = ...
            gonogoPerformance(resp,tt(:,1:2),abort,1,7,[],adjust);
        
        % add to results structure
        s(i).rate = rate;
        s(i).dp = dp;
        s(i).pc = pc;
        s(i).nresp = nresp;
        s(i).ntrials = ntrials;
        
        
        %% psych specific stuff
        if s(i).sesscode == 2
            
            % recompute pc, getting mean first
            s(i).rate = nanmean(s(i).rate,2)';
            s(i).dp = norminv(s(i).rate(2:end)) - norminv(s(i).rate(1));
            s(i).pc = normcdf(s(i).dp / sqrt(2));
         
            
            % fit the psychometric data
            if length(s(i).vols) == 1
                x = p.attenuation;
            else
                x = s(i).vols;
            end
            n = nansum(s(i).ntrials(2:end,:),2);
            if pcFlag
                y = s(i).pc;
            else
                y = s(i).rate(2:end);
            end
            [params,mdl,thresh,sense,fmcon,minfun,thresh75] = ...
                fitLogGrid(x,y,[],n,[],.75);
            
            s(i).prms = params;
            s(i).thresh = thresh;
            s(i).thresh75 = thresh75;
            s(i).sense = sense;
            s(i).maxslope = max(diff(y)./diff(x));
            s(i).fmconFitError = fmcon.minErr;
            s(i).minfunFitError = minfun.fval;
            s(i).minfunParams = minfun.params;
            
        else
            s(i).prms = nan(4,1);
            s(i).thresh = nan;
            s(i).thresh75 = nan;
            s(i).sense = nan;
            s(i).maxslope = nan;
            s(i).fmconFitError = nan;
            s(i).minfunFitError = nan;
            s(i).minfunParams = nan(4,1);
            
        end
        toc;
        
    end
    
    if contains(mouseList{m},'CA123')
    end
    
    tt = struct2table(s);
    if ~iscell(tt.offs)
        tt.offs = num2cell(tt.offs,2);
    end
        if ~iscell(tt.offP)
        tt.offP = num2cell(tt.offP,2);
    end
    
    dat = [dat;tt];
    
end

save(resFile,'dat','mouseList','types');
