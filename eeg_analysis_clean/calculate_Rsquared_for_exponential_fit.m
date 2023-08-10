function [Rsquared] = calculate_Rsquared_for_exponential_fit(params, data, t)
% this function computes the residuals for the exponential fit which is the
% sum of errers data - exp_model
% params is a vector with
%  amplitude, tau, offset


Amp = params(1);
tau = params(2);

model = Amp .* (exp(t/tau)); % compute the model

model(isinf(model)) = 0;
RSS = sum((data - model').^2); % compute RSS
TSS = sum((data - mean(data)).^2); %compute TSS

Rsquared = 1-RSS/TSS;

end 