function [d,res] = muscimol_res(dat)

% mouse and saving options
wd = pwd;
mouseList = {'AW142','AW146','AW147'};
resDir = [wd '/_data'];

% condition for each mouse
condLong = {'0.25mg/mL muscimol',...
                'saline',...
                '0.50mg/mL muscimol'};
condShort = {'musc25','saline','musc50'}; 

% cell analysis options
ops.bin = 0.005;
ops.psthWindow = -.02:ops.bin:.150;
ops.psthTime = ops.psthWindow(1:end-1) + ops.bin/2;
ops.baseLine = [-.1 0];
ops.targetWindow = [.01 .04];

% window analysis options (want to not have overlapping windows)
ops.aucWindow = [ops.psthWindow(1) ops.psthWindow(end)];
ops.timeStep = .0125;
ops.timeWindow = ops.timeStep*2;
ops.timeCent = (min(ops.aucWindow) + ops.timeWindow/2):ops.timeStep:...
    (max(ops.aucWindow) - ops.timeWindow/2);

% population analysis options
ops.noiseLevel = 0;
ops.targetLevel = [5 6];
ops.smooth = .02;
ops.cv = 'leave1out';
ops.respM = {'beh','pop','popCrit','cell','cellCrit'};

% plot options
ops.plotVisible = 'off';
ops.cellPlot = []; %fullfile([wd '/_plots']);
    
cd(wd);

d = dat.d;

if ~exist('res','var')
    
    fprintf('Running muscimol recording analysis for %s... ', ...
            dat.ops.mouse);
    tic;
    
    % run cell analysis
    ops.mouse = dat.ops.mouse;
    ops.condLong = dat.ops.condLong;
    ops.condShort = dat.ops.condShort;
    ops.cellPlot = []; %fullfile([wd '/_plots']); % force to not plot
    [res,ops] = muscimolWrapper(d,ops);
    
end     

cd(wd);
    

    

         

     
     
             
             