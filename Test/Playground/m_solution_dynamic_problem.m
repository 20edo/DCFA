% Compute the solution of the dynamics problem, given the assembled model,
% the time span, the initial conditions (both for displacements and
% velocity) and the vector of loads acting on the structure, N dymention of
% the reduced order problem
%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
function model = m_solution_dynamic_problem(model,tspan,initialcond,f,N)
M = model.M; 
K = model.K; 
C = model.C;
y0 = initialcond; 
[VV,D,FLAG]=eigs(K,M,N,'smallestabs');
w = diag(D).^0.5;
V = VV; 
M_red = diag(diag(V'*M*V)); 
K_red = diag(diag(V'*K*V));
C_red = diag(diag(V'*C*V));
f_red = @(t) V'*f(t); 
 

% solution of the problem
[t,y] = ode45(@(t,y) lin_sys(t,y,M_red,C_red,K_red,f_red), tspan, y0);


% Recovery with modes acceleration 
K_inv = K^-1;
mm = diag(M_red); 
temp = diag(1./(mm.*(w.^2)));  

for i=1:length(t)
   % Recovery of all the displacements  
   U_modes(:,i) = V*y(i,1:N)' + ((K_inv)-V*temp*V')*f(i);
%    U_direct(:,i) = V*y(i,1:n)';
end


for i = 1:length(model.b)
   u_beam = model.b(i).A*U_modes; 
   if i==1 
      model.b(i).in(1).d = zeros(6,1); 
      for k = 2:model.b(i).nel+1
            model.b(i).in(k).d = u_beam(1+6*(i-2):6*(i-1),1);
      end
   else
       for k = 1:model.b(i).nel+1
            model.b(i).in(k).d = u_beam(1+6*(i-2):6*(i-1),1);
       end   
   end
end




end
