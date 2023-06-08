function [res,p] = varianceSims(muTlow,sigmaLow,sigmaHigh,f,sigmaNoise,gainControl,silence)

% load parameters
p = load('optimalNonlinearities_VarianceReconstruction.mat');

fprintf('VOLUME %g\n',muTlow);


%--------------- stimulus parameters and generation ----------------------%

%sd values of the background
%sigmaLow  = 1;
%sigmaHigh = 2;

%target distribution parameters
muThigh = muTlow;
%f = .25;
sigmaTlow  = f .* sigmaLow;
sigmaThigh = f .* sigmaHigh;

%index of the number of spikes (max spikes = levelsInd + 1)
levelsInd = 15;

%width of the variance estimation window
nDecT = 12;

%cycle length 
tC = 50;

%number of switches between states
nC = 1e4;

%total number of samples
nT = tC * nC * 2;

% noise scale
if ~exist('sigmaNoise','var') | isempty(sigmaNoise)
    sigmaNoise = 1;
end

% silent amplitude
if ~exist('silence','var') | isempty(silence)
    silence = 1;
end

%stimulus generation
thetaV = repmat([sigmaLow * ones(tC, 1); sigmaHigh * ones(tC, 1)], nC, 1);
thresh  = (sigmaLow+sigmaHigh)./2;
iL = find(thetaV<thresh);
iH = find(thetaV>=thresh);
xV = thetaV .* randn(nT, 1) .* sigmaNoise .* silence;

%target generation
targetV = zeros(nT,1);
targetV(iL) = randn(numel(iL), 1) * sigmaNoise * sigmaTlow  + muTlow;
targetV(iH) = randn(numel(iH), 1) * sigmaNoise *  sigmaThigh + muThigh;

figure(12343); clf; hold on;
histogram(xV);
histogram(targetV);
drawnow;



%------------------------ allocate results -------------------------------%

%responses to background
respBackgroundV = zeros(nT, 1);

%responses to targets
respTargetV = zeros(nT, 1);

%variance estimate
thetaHatV = zeros(nT, 1);

%decoded value of the stimulus
xHatV = zeros(nT, 1);

%encoding parameter history
kV = zeros(nT, 1);
x0V = zeros(nT, 1);

%decoding parameter history
p0V = zeros(nT, 1);
p1V = zeros(nT, 1);


%---------------------------- run simulation -----------------------------%

%find index of the variance estimate
d = (p.rs_noise01.sig - thetaV(1)).^2;
dInd = find(d == min(d));

%encoding parameters
kCurr = p.rs_noise01.k(dInd, levelsInd);
x0Curr = p.rs_noise01.x0(dInd, levelsInd);

%decoding parameters
p0Curr = p.rs_noise01.p0(dInd, levelsInd);
p1Curr = p.rs_noise01.p1(dInd, levelsInd);


if exist('gainControl','var') & ~isempty(gainControl)
    if strcmp(gainControl,'off')
        kCurr = 0.6;
    elseif strcmp(gainControl,'low')
        kCurr = 0.25;
    end
end


for t = 2:nT
    
     if ~mod(t, 100000)
         fprintf('Step %d out of %d\n', t, nT);
     end
    
    %generate the spiking response to background
    respBackgroundV(t) = encodeWithDiscrete(xV(t), kCurr, x0Curr, p.rs_noise01.sigmaNoise, p.rs_noise01.levels(levelsInd));
    
    %generate the spiking response to the target
    respTargetV(t) = encodeWithDiscrete(targetV(t), kCurr, x0Curr, p.rs_noise01.sigmaNoise, p.rs_noise01.levels(levelsInd));
       
    %decode the background
    xHatV(t) = p1Curr * respBackgroundV(t) + p0Curr;
    
    %update the variance estimate
    if t >= nDecT
        thetaHatV(t) = sqrt(var(xHatV(t-nDecT+1:t)));
    else
        thetaHatV(t) = thetaHatV(t - 1);
    end
    
    %update the decoding parameters
    d = (p.rs_noise01.sig - thetaHatV(t)).^2;
    dInd = find(d == min(d));

    %encoding parameters
    kCurr = p.rs_noise01.k(dInd, levelsInd);
    x0Curr = p.rs_noise01.x0(dInd, levelsInd);

    %decoding parameters
    p0Curr = p.rs_noise01.p0(dInd, levelsInd);
    p1Curr = p.rs_noise01.p1(dInd, levelsInd);
    
    % overwrite kCurr if specified
    if exist('gainControl','var') & ~isempty(gainControl)
        if strcmp(gainControl,'off')
            kCurr = 0.6;
        elseif strcmp(gainControl,'low')
            kCurr = 0.25;
        end
    end
    
    %save parameters
    kV(t) = kCurr;
    x0V(t) = x0Curr;
    p0V(t) = p0Curr;
    p1V(t) = p1Curr;
    
end


res.muT_high  = muThigh;
res.muT_low   = muTlow;
res.sigT_low  = sigmaTlow;
res.sigT_high = sigmaThigh;

res.muB = 0;
res.sigB_high = sigmaHigh; 
res.sigB_low = sigmaLow; 

res.k = kV;
res.x0 = x0V;
res.p1 = p1V;
res.p0 = p0V;
res.yB = respBackgroundV;
res.yT = respTargetV;
res.xB = xV;
res.xT = targetV;
res.xhat = xHatV;
res.theta = thetaV;
res.thetahat = thetaHatV;

res.buffer.k  = buffer(kV(tC+1:end), 3 * tC, tC);
res.buffer.x0 = buffer(x0V(tC+1:end), 3 * tC, tC);
res.buffer.p1 = buffer(p1V(tC+1:end), 3 * tC, tC);
res.buffer.p0 = buffer(p0V(tC+1:end), 3 * tC, tC);
res.buffer.yB = buffer(respBackgroundV(tC+1:end), 3 * tC, tC);
res.buffer.yT = buffer(respTargetV(tC+1:end), 3 * tC, tC);
res.buffer.xB = buffer(xV(tC+1:end), 3 * tC, tC);
res.buffer.xT = buffer(targetV(tC+1:end), 3 * tC, tC);
res.buffer.xhat = buffer(xHatV(tC+1:end), 3 * tC, tC);
res.buffer.theta = buffer(thetaV(tC+1:end), 3 * tC, tC);
res.buffer.thetahat = buffer(thetaHatV(tC+1:end), 3 * tC, tC);

%---------- compute entropy of the response and discriminability ---------%


%temporal indices
tInds = 1:3 * tC;

%allocate results
hC = 0:p.rs_noise01.levels(levelsInd)-1;
klD = zeros(length(tInds), 2);

%buffer data
rB = buffer(respBackgroundV(tC+1:end), 3 * tC, tC);
rT = buffer(respTargetV(tC+1:end), 3 * tC, tC);

for t = 1:length(tInds)

    cnB = histc(rB(t, :), hC) + 1;
    cnT = histc(rT(t, :), hC) + 1;
    
    cnB = cnB ./ sum(cnB);
    cnT = cnT ./ sum(cnT);
    
    %KL-divergence or symmetrized KL-divergence
    indsP = cnB > 0 & cnT > 0;
    klD(t, 2) = 0.5 * (sum(cnB(indsP) .* log2(cnB(indsP) ./ cnT(indsP))) + sum(cnT(indsP) .* log2(cnT(indsP) ./ cnB(indsP))));
   
    %overlap between distributions
    klD(t, 1) = 1 - sqrt(cnB * cnT');

end

res.buffer.Dkl = klD;


end