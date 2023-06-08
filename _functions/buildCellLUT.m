clear all; close all;

% file list
dataDir = '~/data/gain_behavior/_cells';
files = dir(fullfile(dataDir,'*.mat'));
files(contains({files.name},'sim') | contains({files.name},'LUT')) ...
    = [];

T = struct;

% load each file
for i = 1:length(files)
    
    fprintf('%d/%d ',i,length(files)); tic;
    
    load(fullfile(files(i).folder,files(i).name));
    
    if isnan(d.waveform)
        waveform = nan(1,63);
    else
        waveform = d.waveform;
    end
    wave = normalize(waveform,'range',[-1 1]);
    waveformNorm = wave - mean(wave(1));

    % build the table
    T(i).mouse = d.mouse;
    T(i).session = d.session;
    T(i).cellcnt = d.cellcnt;
    T(i).cellnum = d.cellnum;
    T(i).celltype = d.cellType;
    T(i).cellID = d.cellID;
    T(i).FR = length(d.spikes) ./ (d.spikes(end) - d.spikes(1));
    T(i).NR = nanmean(d.NR,1);
    T(i).fwhm = d.FWHM;
    T(i).trough_peak = d.trough_peak;
    T(i).peak_inflect = d.peak_inflect;
    T(i).snr = d.snr;
    T(i).ISIcontamination = d.ISIcontamination;
    T(i).cellFile = fullfile(files(i).folder,files(i).name);
    T(i).stimFile = d.stimFile;
    T(i).behaviorFile = d.behaviorFile;
    T(i).neuralFile = d.neuralFile;
    T(i).sd = stimInfo.sd;
    T(i).contrast = stimInfo.sd(2) > stimInfo.sd(1);
    T(i).waveform = waveform * .195;  % convert to volts
    T(i).waveformNorm = waveformNorm;
    
    toc;
    
    
end

t = struct2table(T);

%% group neurons by fast or regular spiking
wavecorr = corr(nanmean(t.waveform)',t.waveform');
cutoff = .5;

% cluster all units
X = [t.trough_peak(wavecorr>=cutoff), ...
     t.fwhm(wavecorr>=cutoff), ...
     t.FR(wavecorr>=cutoff)];
dim2use = 1;
k = 2;
gm = fitgmdist(X(:,dim2use),k,'CovarianceType','full','SharedCovariance',false);
clusterx = cluster(gm,X(:,dim2use));
dists = pdist([zeros(1,length(dim2use)); gm.mu]);

% order by which is closer to origin (FS are closer)
if dists(2) > dists(1)
    remap = [2 1];
    clusterx = remap(clusterx)';
end 
kcols = [0 0 0; 1 0 1];
kcv = kcols(clusterx,:);

% place cluster index into original index
t.spiketype = zeros(size(t,1),1);
t.spiketype(wavecorr>=cutoff) = clusterx;

% plot histogram of correlations
clf; 
subplot(2,3,1); hold on;
edges = linspace(-1,1,21);
histogram(wavecorr(wavecorr>=cutoff),edges,'facecolor',[0 0 0]);
histogram(wavecorr(wavecorr<cutoff),edges,'facecolor',[.5 .5 .5]);
plot([.5 .5],ylim,'r--','linewidth',1);
xlabel('Correlation with mean waveform');
ylabel('Count'); plotPrefs; axis tight;

% plot good and bad waveforms
subplot(2,3,2); hold on;
plot(t.waveform(wavecorr>=cutoff,:)','color',[0 0 0]);
plot(t.waveform(wavecorr<cutoff,:)','color',[.5 .5 .5]);
xlabel('Time (bins)'); ylabel('Voltage (uV)');

% scatter plot of all units + FR
subplot(2,3,3); hold on;
scatter(t.trough_peak(wavecorr>=cutoff)*1000,...
        t.fwhm(wavecorr>=cutoff)*1000,10,t.FR(wavecorr>=cutoff),'filled')
ch = colorbar; ylabel(ch,'spks/s'); xlim([0 1]); ylim([0 .5]);
xlabel('Trough-Peak (ms)'); ylabel('FWHM (ms)');

% labeled scatter
subplot(2,3,6); hold on;
%scatter3(X(:,1),X(:,2),X(:,3),10,kcv,'filled')
scatter(X(:,1)*1000,X(:,2)*1000,10,kcv,'filled');
ch = colorbar; xlim([0 1]); ylim([0 .5]);
xlabel('Trough-Peak (ms)'); ylabel('FWHM (ms)');

% labelled waveforms (normalized)
subplot(2,3,5); hold on;
ph = patchErrorBars(1:63,t.waveform(t.spiketype == 1,:),...
                    kcols(1,:),'std');
ph = patchErrorBars(1:63,t.waveform(t.spiketype == 2,:),...
               kcols(2,:),'std');

% save out the table
save(fullfile(dataDir,'..','cellLUT.mat'),'t');








    
    











