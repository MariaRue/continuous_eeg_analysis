function model = calculate_exponential_model_for_estimated_parameters(params,t)
Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model
end