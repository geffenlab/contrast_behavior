function session = matchBehaviorToSpikes(mouseList,behaviorBase,spikeBase,onsetDiff)

%% function session = matchBehaviorToSpikes(mouseList,behaviorBase,spikeBase,onsetDiff)
%
% This function matches neural recording sessions to their
% corresponding behavioral sessions by date and time. 
%  INPUTS
%   mouseList:
%   behaviorPath:  path to a mouse behavior folder, containing all sessions
%   spikePath:     path to the same mouses spike folder
%   onsetDiff:    threshold in minutes for session onset separation.
%                 if both start within this time, the first session was
%                 likely aborted, and will not be considered (default 30)

if nargin < 3
    onsetDiff = 30;
end

fprintf('Matching neural and behavioral recordings:\n');

% setup behavior types
types = {'noresponse','training','psychometric','offset','opto-psychometric','threshold-opto','switch'};

% for each mouse
nSessions = 0;
for m = 1:length(mouseList)
    
    mouse = mouseList{m};
    behaviorDir = fullfile(behaviorBase,'_data',mouse);
    spikeDir = fullfile(spikeBase,mouse);
    
    fprintf('%s\n',mouse);
    
    % find behavioral files
    [fileList fileInd] = indexDataFiles(behaviorDir); 
    ind = fileInd(:,5) == 1 & fileInd(:,2) > 1;
    fileInfo = fileInd(ind,:);
    files = fileList(ind);
        
    if exist(spikeDir,'dir') & ~isempty(dir(spikeDir))

        for i = 1:size(fileInfo,1)
            
            % format the date
            dt = datetime(num2str(fileInfo(i,6)),'InputFormat','yyMMddHHmm');
            
            % get a list of spike sessions directories
            d = dir(spikeDir); d(1:2) = []; d = d([d.isdir]);
            
            % discard directories that aren't recording length
            d = d(cellfun(@length,{d.name})==19);
            
            for j = 1:length(d)
                dt_spike(j) = datetime(d(j).name,'InputFormat',...
                                       'yyyy-MM-dd_HH-mm-ss');
                
            end
            
            % find onset differences below specified number
            ind0 = abs(dt_spike - dt) < minutes(onsetDiff);
            
            if sum(ind0) > 1
                fprintf(['\tMore than one spike sessions was started within ' ...
                         '30 minutes of this behavioral session!\n']);
                
                % check for sorted data
                fI = find(ind0);
                ind = ind0;
                cnt = 0;
                for ii = find(ind)
                    
                    df = dir(fullfile(d(ii).folder,d(ii).name,...
                                      'experiment1','recording1', ...
                                      'continuous','Rhythm_FPGA-100.0'));
                    if ~any(contains({df.name},'cluster_group'))
                        ind(ii) = 0;
                    else
                        cnt = cnt + 1;
                        spkFile = contains({df.name},'spike_times.npy');
                        spikeSize(cnt) = df(spkFile).bytes;
                        
                    end
                    
                end
                
                % if there are still multiple options, choose the largest
                fI = find(ind);
                if length(fI) > 1
                    mi = spikeSize == max(spikeSize);
                    ind(fI(~mi)) = 0;
                    
                end
                
            else 
                ind = ind0;
                
            end

            if sum(ind) == 1
                nSessions = nSessions + 1;
                session(nSessions).mouse = mouse;
                session(nSessions).behavior = files{i};
                session(nSessions).neural = fullfile(d(ind).folder,d(ind).name);
                session(nSessions).session = datestr(dt_spike(ind),...
                                                     'yymmdd');
                session(nSessions).type = types{fileInfo(i,2)+1};
                session(nSessions).sessionInfo = fileInfo(i,:);
                
            elseif sum(ind) > 1
                fprintf(['More than one spike sessions was started within ' ...
                         '30 minutes of this behavioral session!\n']);
                
                % check if each folder has 
                keyboard
                
            end
            
        end
        
    end

end

fprintf('\n');