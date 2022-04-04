function [residual_reg] = calculate_residuals_for_exponential_fit(params, data, t)
% this function computes the residuals for the exponential fit which is the
% sum of errers data - exp_model
% params is a vector with
%  amplitude, tau, offset


Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model

model(isinf(model)) = 0;
residual = sum((data - model').^2); % compute the error
residual_reg = residual + sum(params.^2) .* 0.01;
end 