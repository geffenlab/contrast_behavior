clear all; close all;
%addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('./_functions/'));

% paths
root = pwd;
dataDir = './_data'; % location of data files
addpath(genpath(dataDir));

%% figure 2: GC-GLM
% to run this, clone https://github.com/chris-angeloni/contrast_glm
cd(root);
resPath = fullfile(root,'_data','_glm','_res');
cd ../contrast_glm/ % replace with path to contrast_glm repo here
run_gcglm; cd(root); clear ops res r;


%% behavioral analysis setup
% load all spike and session data
fprintf('Loading spike data... '); tic;
spikeData = load(fullfile(dataDir,'spikeData.mat')); toc;
fprintf('Loading session data... '); tic;
sessionData = load(fullfile(dataDir,'sessionData.mat')); toc;

% set the seed
ops.seed = 1989;
rng(ops.seed);

% waveform analysis
res_wave = run_wave(spikeData);
%f1 = figure(1); clf; plot_waveforms(res_wave);

% cell inclusion criteria
%  1) Spike rate > 1Hz
%  2) good waveform shapes (not inverted, good pvalue, not super wide)
wave_include = ones(size(spikeData.cellinfo,1),1);
wave_include(~isnan(vertcat(spikeData.waveform.FWHM))) = res_wave.include;
include = [spikeData.cellinfo{:,8}]'>1 & wave_include;
ops.include = include;


%% Figure 3: run and plot behavior
% (requires behavior folder: ~/gits/gain-gonogo/')
stats_beh = run_behavior;


%% Figure 1: run and plot the normative model simulations
% (this takes 20-30 minutes the first time, several minutes after
%  saving the initial simulatuion results)
stats_norm = run_normative_model;

% plot percent change in behavior versus model
fh = plot_behavior_model_comp(stats_beh,stats_norm,1);
saveFigPDF(fh,[600 200],'./_plots/_behavior_model_comp.pdf');



%% Figure 4: run and plot muscimol results
% (requires behavior folder: ~/gits/gain-gonogo/')
% (this will take roughly 20 minutes first run, then ~5 min rerunning)
stats_muscimol = run_muscimol;


%% Figure 5: psych
% (this will take quite a while the first time, several hours)

% options
ops.resDir = './_data';
ops.bin = .005;
ops.baseline = [-.1 0];
ops.target = [0 .1];
ops.response = [.1 1];
ops.edges = ops.target(1)+ops.baseline(1) : ...
    ops.bin : ...
    ops.target(2)+ops.response(2)-ops.baseline(1);
ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
ops.targetLevel = [5 6];
ops.noiseLevel = [0];
ops.smooth = 2;
ops.timeInd = 1;
ops.sig_neurons = true;
ops.resFile = '_res_psych.mat';

% run
[res_psych,r_psych] = run_psych(spikeData,sessionData,ops);

% plot summary
stats_psych = plot_psych_summaries(r_psych,ops);

% plot example neuron
f41 = figure(41); clf;
%res_psych.single_cell(find(contains({res_psych.single_cell.cellID},'CA118_2007070944'))).cellID
cellID = 'CA118_2007070944_09_040_mu';
plot_single_cell_psych(cellID,res_psych,spikeData,sessionData,ops);
saveFigPDF(f41,[300 475],'./_plots/_psych_ex_neuron_loc.pdf');

% plot example sessions
sz = [500 520];
f42 = figure(42); clf; set(f42,'Position',[0 0 sz]);
sessID = 2007070944; % ca118 - low contrast
plot_session_summary_psych_fig(sessID,res_psych,spikeData,sessionData,ops,...
                               r_psych,'critp_adj');
saveFigPDF(f42,sz,'./_plots/_psych_ex_session_loc.pdf');

f43 = figure(43); clf; set(f43,'Position',[0 0 sz]);
sessID = 2006110945; % ca118 - high contrast
plot_session_summary_psych_fig(sessID,res_psych,spikeData,sessionData,ops,...
                               r_psych,'critp_adj');
saveFigPDF(f43,sz,'./_plots/_psych_ex_session_hic.pdf');




%% Figure 5: offset

% options
ops.resDir = './_data';
ops.bin = .005;
ops.baseline = [-.1 0];
ops.target = [0 .1];
ops.response = [.1 1];
ops.edges = ops.target(1)+ops.baseline(1) : ...
    ops.bin : ...
    ops.target(2)+ops.response(2)-ops.baseline(1);
ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
ops.targetLevel = [2];
ops.noiseLevel = [0];
ops.smooth = 2;
ops.timeInd = 1;
ops.include = include;
ops.sig_neurons = false;
ops.resFile = '_res_offset.mat';

% run
[res_off,r_off] = run_offset(spikeData,sessionData,ops);

% plot summary
stats_off = plot_offset_summaries(r_off);

% example neuronmake
plot_single_cell_offset(...
    spikeData.cellinfo{26,7},res_off,spikeData,sessionData,ops);



%% figure 6: LN model
% (fitting each cell takes some time... ~2-3 days for the initial run)

% model options
ops.bin = .025;
ops.fs = 1/ops.bin;
ops.w = 0.3;
ops.pad = ops.w;
ops.f = []; % frequencies depend on the data
ops.t = (-ops.w:1/ops.fs:0);
ops.nbins = 50;
ops.weight = true;
ops.sigma = .01;
ops.includeZero = false;
ops.modelType = 'exponential';
ops.kfolds = 10;
ops.models = {'sta','NRC','glm'};
ops.modelStrings = {'Spike-Triggered Average',...
                    'Normalized Reverse Correlation',...
                    'Regularized GLM'};
ops.targetExclude = .05;
ops.trialCutoff = 4.5 * ops.fs;
ops.resDir = './_data';
ops.resFile = '_res_lnmodel.mat';
ops.fig_visible = 'off';

% run
[res_ln] = run_lnmodel(spikeData,sessionData,ops);

% plots and stats
stats_ln = plot_lnmodel_summaries(res_ln,res_psych,r_psych,ops);
