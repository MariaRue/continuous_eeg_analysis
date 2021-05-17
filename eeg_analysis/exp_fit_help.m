lags = 0:100;

tau = 100;
t_offset = 34.2;

model(lags<t_offset) = 0;
model(lags>=t_offset) = exp(-(lags(lags>=t_offset)-t_offset)/tau);

model = model/sum(model);

plot(model);hold on;


t_offset = 35.5;

model(lags<t_offset) = 0;
model(lags>=t_offset) = exp(-(lags(lags>=t_offset)-t_offset)/tau);

model = model/sum(model);

plot(model);

%model = A*model;
