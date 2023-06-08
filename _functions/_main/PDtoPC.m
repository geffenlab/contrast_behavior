function pc = PDtoPC(hr,fa)

%% function pc = PDtoPC(hr,fa)
%
% This function converts performance in terms of p(detect) or
% p(respond) to percent correct (pc), by normalizing d'.

hr(hr==1) = .999;
hr(hr==0) = .001;
fa(fa==1) = .999;
fa(fa==0) = .001;

pc = normcdf( (norminv(hr) - norminv(fa)) ./ sqrt(2) );