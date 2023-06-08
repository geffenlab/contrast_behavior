function [pred, crit, rule] = crit_classifier(spikes,index)

%% function [pred, crit, rule] = crit_classifier(spikes,index,tt)
%
% Runs a simple classifier that first finds the criterion that best
% discriminates between spikes where index == 0 and spikes where
% index ~= 0 (does this twice, checking for a greater than or less
% than decision rule)

% set up criterion
mn = min(spikes);
if all(spikes > 0) & mn ~= 0
    mn = 0;
end
mx = max(spikes);
if all(spikes < 1) & mx == 0
    mx = 1;
end
crits = repmat(linspace(mn,mx,100),1,2);

% find max performance crit
for i = 1:length(crits)
    if i <= 100
        % say detect if spikes < crit
        pred = spikes < crits(i);
    else
        % say detect if spikes > crit
        pred = spikes > crits(i);
    end
    critpc(i) = mean(pred == (index(:,1) ~= 0),'omitnan');
end

% find max performance
[~,mi] = max(critpc);
crit = crits(mi);

% determine decision rule
if mi <= 100
    rule = @(spikes,crit) spikes < crit;
else
    rule = @(spikes,crit) spikes > crit;
end

% using that crit, generate spiking prediction of target vs noise
pred = rule(spikes,crit);
