% understand what is meant by 'OFFSET in the Murry paper' 

x = 0.1:0.001:0.51; 

B = 0; 
tau = 0.3; 
A = 1;


y = A .* (exp(x/tau)+B);

plot(y,'b')

hold on 
A = 1; 
B = 0; 
y = A .* (exp((-x/tau))+B);
plot(y,'k')

hold on 
A = 3; 
B = 0.4; 


y = A .* (exp(-x/tau)+B);
plot(y,'g')
