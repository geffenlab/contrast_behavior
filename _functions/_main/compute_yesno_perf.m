function [r,ny,nt] = compute_yesno_perf(response,index,conditions,adjust)

%% function [r,ny,nt] = compute_yesno_perf(response,index,conditions,adjust)
%
% This function takes responses, a trial index, and desired unique
% conditions and computes for each condition the number of trials, 
% the number of yes responses, and the response rate
% 
% optionally adjusts the data using the log-linear rule (Hautus,
% 1995)

% if no conditions provided, use unique rows in index
if ~exist('conditions','var') | isempty(conditions)
    uI = unique(index,'row');
else
    uI = unique(conditions,'row');
end

% optional adjustments
if ~exist('adjust','var') | isempty(adjust)
    adjust = false;
end

if adjust
    % log-linear adjustment
    adj = [.5 1];
else
    % no adjustment
    adj = [0 0];
end

% compute response/trial counts and rates
ny = nan(1,length(uI));
nt = nan(1,length(uI));
r  = nan(1,length(uI));
for i = 1:length(uI)
    
    I = all(index == uI(i,:),2);
    
    if sum(I) > 0
        ny(i) = sum(response(I));
        nt(i) = sum(I);
        r(i) = (ny(i) + adj(1)) / (nt(i) + adj(2));
    end
end


