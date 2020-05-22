clear all 
cd .. 
Clamped_eigenshapes

C = 0.2*K; 

f = @(t) [zeros(size(M,1)-6,1); 1e6*(t>0); zeros(5,1)];
y0 = zeros(2*n,1);
M_red = diag(diag(V'*M*V)); 
K_red = diag(diag(V'*K*V));
C_red = diag(diag(V'*C*V));
f_red = @(t) V'*f(t); 
cd Test
cd Playground
% time interval
tspan = [-1,1000]; 
% initial conditions
n = size(M_red,1); 



[t,y] = ode45(@(t,y) lin_sys(t,y,M_red,C_red,K_red,f_red), tspan, y0);

% U = @(t) V*y(:,1:n)' + ((K^-1)-V*(M_red^-1*(diag(w.^2))^-1)*V')*f(t);

K_inv = K^-1;
mm = diag(M_red); 
temp = diag(1./(mm.*(w.^2)));  
for i=1:length(t)
   % vedere lo smorzamento
   U_modes(i,:) = (V*y(i,1:n)' + ((K_inv)-V*temp*V')*f(i))';
%    figure(1)
%    plot(6:6:size(M,1),U(i,6:6:end))
%    pause(0.05)
%    U_direct(i,:) = (V*y(i,1:n)')';
end

for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*U_modes; 
    for k = 1:model.b(i).nel+1
       model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),end);
%        if we would like to reconstruct the solution over time
%        model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end








