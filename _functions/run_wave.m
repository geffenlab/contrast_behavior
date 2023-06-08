function res = run_wave(spikeData)

% concatenate all waveforms and normalize
waveforms = vertcat(spikeData.waveform.peakwaveform);
nwaveforms = normalize(waveforms,2);

% time vector
t0 = (1:size(waveforms,2)) / 30000;

% concatenate waveform metrics, remove nan
trough_peak = vertcat(spikeData.waveform.trough_peak);
trough_peak(isnan(trough_peak)) = [];
FWHM = vertcat(spikeData.waveform.FWHM);
FWHM(isnan(FWHM)) = [];
nwaveforms(all(isnan(nwaveforms),2),:) = [];

% correlate each waveform to the mean
[r,p] = corr(nwaveforms',mean(nwaveforms,1)');

% cut off low correlated waveforms
rcut = r > 0;

% sort by correlation
[rsort,sorti] = sort(r,'ascend');
psort = p(sorti);
[~,r0] = min(abs(rsort));

% corrected p value, get the index
pcorr = .05 / length(r);
pcut = p < pcorr;
dp = abs(psort-pcorr);
dppos = dp(rsort>0);
dpmin = min(dppos);
pcuti = find(abs(psort-pcorr) == dpmin);

% remove FWHM outliers (prolly noise units)
outliers = isoutlier(FWHM);

% inclusion index
include = pcut & rcut & ~outliers;

% cluster the remaining waveforms
X = [trough_peak(include), FWHM(include)]; k = 2;
gm = fitgmdist(X,k,'CovarianceType','full','SharedCovariance',false);
clusterx = cluster(gm,X);
dists = pdist([0 0; gm.mu]);


% order by which is closer to origin (FS are closer)
if dists(2) > dists(1)
    remap = [2 1];
    clusterx = remap(clusterx)';
end 


% results
res.waveforms = waveforms;
res.nwaveforms = nwaveforms;
res.r = r;
res.p = p;
res.r0 = r0;
res.pcuti = pcuti;
res.rcut = rcut;
res.pcut = pcut;
res.outliers = outliers;
res.include = include;
res.good_waveforms = include;
res.type = clusterx;
res.X = X;
res.t0 = t0;