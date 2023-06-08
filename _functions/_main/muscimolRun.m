clear all
close all

warning off

%% directory stuff
spikeBase = '~/data/kilosort';
wd = '~/chris-lab/projects/contrast_behavior';
addpath(genpath('~/gits/npy-matlab'));
addpath(genpath('~/chris-lab/code_general'));

% mouse and saving options
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
ops.cellPlot = fullfile([wd '/_plots']);

% run options
dryRun = false;

% for each mouse
for i = 1:length(mouseList)
    
    fn = fullfile(resDir,...
                  sprintf('%s_%s.mat',mouseList{i},condShort{i}));
    
    % build blocked data if needed
    if ~exist(fn,'file')
        r = dir(fullfile(spikeBase,mouseList{i},'**','Events.nev'));
        root = r.folder;
        load('blockTemplates.mat');
        [d.block, d.cellInfo, d.labels, d.spikes, d.events, d.fs] = ...
            splitSpikesByTemplates_MSE(root,template,true,'NLX');
        save(fn,'d');
    else
        load(fn);
    end
    
    if ~exist('res','var')
        cd(wd);
        
        % run cell analysis
        ops.mouse = mouseList{i};
        ops.condLong = condLong{i};
        ops.condShort = condShort{i};
        [res,ops] = muscimolWrapper(d,ops);
    
        save(fn,'res','ops','-append');
        clear res
        
    else
        % run population analysis
        %[res,ops] = muscimolPop(res,ops);
        
    end
    
end

cd(wd);
    
% plot style
ls = {'--','-','-.'};
ms = 4;
mstyle = {'^','o','x'};
colors = {'b','r'};
snrs = d.block(1).stimInfo.p.targetDBShift;
xtick = [-5 snrs];
xtickl = num2str([-inf snrs]');

f1 = figure(1); clf;
ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,2);
ax(3) = subplot(2,2,3);
ax(4) = subplot(2,2,4);

% compare effects over different sessions
order = [2 1]; clear hh;
for i = 1:length(order)
    
    I = order(i);
    
    % load data
     fn = fullfile(resDir,...
                  sprintf('%s_%s.mat',mouseList{I},condShort{I}));
     load(fn)
     
     % for each cell
     for j = 1:length(res)
         % extract post session
         targMean(j,:,:) = mean(res(j).condTarget{2});
         targAUC(j,:,:) = res(j).auc;
     end
     
     % plot each contrast
     clear x y;
     for j = 1:2
         axes(ax(j)); hold on;
         
         % fit
         x = xtick(2:end);
         y = mean(squeeze(targMean(:,j,2:end)));
         [prms,mdl,thresh] = fitLogGrid(x,y);
         
         % plot
         h = errorBars(snrs,squeeze(targMean(:,j,2:end)),colors{j});
         h.Marker = mstyle{I}; h.LineStyle = ls{I};
         hn = errorBars(-5,squeeze(targMean(:,j,1)),colors{j});
         hn.Marker = mstyle{I};
         
         axes(ax(j+2)); hold on;
         
         % fit
         x = xtick(2:end);
         y = mean(squeeze(targAUC(:,j,:)));
         [prms,mdl,thresh] = fitLogGrid(x,y);
         
         % plot
         h = errorBars(snrs,squeeze(targAUC(:,j,:)),colors{j});
         h.Marker = mstyle{I}; h.LineStyle = ls{I};

     end
     
     hh(I) = h;
end
legend(hh,condLong(1:length(order)))

yl = get(ax,'ylim');
ym = [min([yl{:}]) max([yl{:}])];
for i = 1:2
    set(ax(i),'ylim',ym,'xtick',xtick,'xticklabels',xtickl);
    axes(ax(i)); xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
    plotPrefs;
end
for i = 3:4
    set(ax(i),'ylim',ym,'xtick',xtick,'xticklabels',xtickl);
    axes(ax(i)); xlabel('Target Volume (dB SNR)'); ylabel('spks/s');
    plotPrefs;
end
%saveFigPDF(f1,[500 200],'./_plots/_figures/sfig_muscimol_dose_actx.pdf');
    

         

     
     
             
             