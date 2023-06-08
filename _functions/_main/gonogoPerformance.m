function [rate, dp, pc, nresp, ntrials] = gonogoPerformance(resp,tt,abort,start,crit,dim,adjust)

%% function [rate, dp, pc, nresp, ntrials] = gonogoPerformance(resp,tt,abort,start,crit,dim,adjust)
%
% this function computes performance in a go/no-go task, where
% there are either nogo trials, where the correct response is 0, or
% go trials where the correct response is 1. 
%
% INPUTS: 
%  resp = response vector, [0] for no response [1] for response
%  tt = trial vector, which tells the condition for each trial.
%       (tt[:,1] is interpreted as the target condition, where 0
%        is interpreted as no-go, and above 0 is go)
%  abort = abort vector, doesn't count these trials
%  start = exludes 1:start trials
%  crit = excludes find(sum(resp) - cumsum(resp) == crit):end trials
%  dim = optionally average performance over one dimension 
%        (1 to average over volume, 2 to average over time)
%
% OUTPUTS:
%  rate = hit rate (or FA if 0)
%  dp = dprime relative to 0 trials (length will be -1 rate)
%  pc = percent correct relative to 0 trials
%  nresp = response frequency per condition
%  ntrials = ntrials per condition

if ~exist('abort','var') | isempty(abort)
    abort = zeros(size(resp));
end

% first normalize inputs
if size(resp,2) > size(resp,1)
    resp = resp';
end
if size(tt,2) > size(tt,1)
    tt = tt';
end
if size(abort,2) > size(abort,1)
    abort = abort';
end

if size(tt,2) > 2
    error(['gonogoPerformance: This function only takes 2 columns ' ...
           'for trial index!!']);
end

% trim early aborted trials
mn = min([size(resp,1) size(tt,1) size(abort,1)]);
resp = resp(1:mn);
tt = tt(1:mn,:);
abort = abort(1:mn);

% find good trials where mouse is performing
if ~exist('start','var') | isempty(start)
    start = 1;
end

endv = mn;
if exist('crit','var') & ~isempty(crit)
    % count back response to crit and include that trial
    endv = find((sum(resp) - cumsum(resp)) == crit);
    
end

% mask trials
goodInd = zeros(mn,1);
goodInd(start:endv) = 1;
goodInd = logical(goodInd);

% compute performance for each unique value in tt
u1 = unique(tt(:,1));
u2 = unique(tt(:,2));

for i = 1:length(u1)
    for j = 1:length(u2)
        
        I = all(tt == [u1(i) u2(j)],2);
        nresp(i,j) = sum(resp(I));
        ntrials(i,j) = sum(I);
        
    end
end

% optional adjustments
if ~exist('adjust','var') | isempty(adjust)
    adjust = '';
end

if contains(adjust,'log')
    % log-linear adjustment
    rate = (nresp+.5) ./ (ntrials+1);

elseif contains(adjust,'2N')
    % 2N adjustment
    rate = nresp ./ ntrials;
    rate(rate == 0) = 1 ./ (2*ntrials(rate==0));
    rate(rate == 1) = 1 - (1 ./ (2*ntrials(rate==1)));
else   
    % arbitrary rate adjustment
    rate = nresp ./ ntrials;
    rate(rate == 0) = .001;
    rate(rate == 1) = .999;
end

% average if specified
if exist('dim','var') & ~isempty(dim)
    rate = nanmean(rate,dim);
    nresp = sum(nresp,dim);
    ntrials = sum(ntrials,dim);
end

% d-prime and percent correct
dp = norminv(rate(2:end,:)) - norminv(rate(1,:));
pc = normcdf(dp ./ sqrt(2));