clear all
close all
warning off


%% directory stuff
behaviorBase = '~/gits/gain-gonogo';
spikeBase = '~/data/kilosort';
wd = '~/chris-lab/projects/gain_behavior';
addpath(genpath(behaviorBase));
addpath(genpath('~/gits/npy-matlab'));
addpath(genpath('~/chris-lab/code_general'));
addpath(genpath([wd filesep]));
rmpath([wd '/!archive']);

% mouse and saving options
mouseList = {'CA102','CA104','CA105','CA106','CA107',...
            'CA118','CA119','CA121','CA122'};
datDir = [wd '/_data'];
resDir = [wd '/_data/GLM'];
cellPlot = fullfile([wd '/_plots']);


%% model options
ln.bin = .025;
ln.fs = 1/ln.bin;
ln.w = 0.3;
ln.f = []; % frequencies depend on the data
ln.t = (-ln.w:1/ln.fs:0);
ln.nbins = 50;
ln.weight = true;
ln.sigma = .01;
ln.includeZeros = true;
ln.modelType = 'exponential';
ln.kfolds = 10;
ln.models = {'sta','nrc','glm'};
ln.modelStrings = {'Spike-Triggered Average',...
                   'Normalized Reverse Correlation',...
                   'Regularized GLM'};


%% run options
dryRun = false;
debugMode = false;
plotVisible = 'off';
plotSize = [1500 400];
resAppend = '_GLM';

if debugMode
    fprintf('\n***\nWARNING: RUNNING IN DEBUG MODE!!\n***\n');
    resAppend = '_debug';
end


%% build a session list
session = matchBehaviorToSpikes(mouseList,behaviorBase,spikeBase,30);
%session = session(contains({session.mouse},'CA121') & contains({session.session},'200413'));




%% session loop
for i = 1:length(session)
    
    fprintf('MOUSE %s - SESSION %s\n',session(i).mouse,session(i).session);
    fnout = makeSessionData(session(i).mouse,session(i),datDir);
    
    % skip files with no data
    if isempty(fnout)
       continue;
    end
    
    if ~dryRun
        
        % results file
        [~,fn] = fileparts(fnout);
        resFile = fullfile(resDir,[fn resAppend '.mat']);
                
        res = struct([]);
        if ~exist(resFile,'file')
            
            % start figure
            f1 = figure(1); set(f1,'visible',plotVisible);
            saveFigPDF(f1,plotSize);
            
            % load the raw data
            load(fnout);
            
            % frequency
            ln.f = data.behavior.stimInfo.freqs;
            
            % cell loop
            for c = 1:data.params.nCells
                
                fprintf('\tCell %d/%d\n',c,data.params.nCells);
                
                
                %% fit the LN model
                [res,data,ops,exitflag] = lnmodelcv(data,res,c,cellPlot,ln,debugMode);
                
                
                % save plot
                if exist('cellPlot','var') & ~isempty(cellPlot) & ~exitflag
                    tStr = sprintf('%s-cell%d-%s%s',data.params.plot.titleString,...
                                   c,data.cellInfo{c,6},resAppend);
                    fn = sprintf('./_plots/%s.pdf',tStr);
                    saveFigPDF(f1,plotSize,fn)
                    clf(f1);
                end
                
            end
            
            % save res and data
            save(resFile,'res','data','ops');
            
        end
        
    end
    
    fprintf('-----------------\n\n');
    
end

%gainRes

