function [errr,y_pred] = gain_fun(tau,y,X,contrast,coeffs,obj,ops)

% linear drive
linear_drive = coeffs(1) + X * coeffs(2:end);

% true tau
tau_true = obj.tau;

% previous method using model parameters:
% x(x==0) = obj.x0;
% x0 = convSTRF(x'-obj.operating_point,fliplr(obj.beta))';
% linear_drive = log(obj.base_rate) + x0;

% simulate gain control with new tau
obj.tau = tau;
[drive_with_gain,g] = gain_scaling(obj, linear_drive, contrast, ops);

% exponentiate for poisson rate, generate poisson spikes
y_pred = exp(drive_with_gain);
%y_pred = poissrnd(l);

% error (only around the transition)
errr = norm(y(ops.transition)-y_pred(ops.transition));