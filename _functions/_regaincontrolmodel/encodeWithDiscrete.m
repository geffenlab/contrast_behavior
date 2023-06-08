function [rV,rV0] = encodeWithDiscrete(x, k, x0, sigmaNoise, nLevels)

%response vector
rV = 1 ./ (1 + exp(-k * (x - x0) ) );
rV0 = rV;

%Gaussian noise
rV = rV + randn(size(rV)) * sigmaNoise;
rV(rV<0)=0;
rV(rV>1)=1;

bins = linspace(0,1,nLevels+1);
spikes = (0:(nLevels-1))';

[~,~,binCounts] = histcounts(rV,bins);
rV = spikes(binCounts);
    