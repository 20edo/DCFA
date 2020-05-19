clear all 
cd .. 
Clamped_eigenshapes

C = 0.1*K; 

f = @(t) [zeros(size(M,1)-1,1); 1e6*sin(t/(2*pi))];
M_red = diag(diag(V'*M*V)); 
K_red = diag(diag(V'*K*V));
C_red = diag(diag(V'*C*V));
f_red = @(t) V'*f(t); 
cd Test
cd Playground
% time interval
tspan = [0,10]; 
% initial conditions
n = size(M_red,1); 
y0 = zeros(2*n,1);


[t,y] = ode45(@(t,y) lin_sys(t,y,M_red,C_red,K_red,f_red), tspan, y0);

% U = @(t) V*y(:,1:n)' + ((K^-1)-V*(M_red^-1*(diag(w.^2))^-1)*V')*f(t);


for i=1:length(t)
   U(i,:) = (V*y(i,1:n)' + (real(K^-1)-V*(real(M_red^-1)*(diag(w.^2))^-1)*V')*f(i))';
   figure(1) 
   plot(6:6:size(M,1),U(i,6:6:end))
   pause(0.05)
end