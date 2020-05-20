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
C = K;
y0 = initialcond; 
[VV,D,FLAG]=eigs(K,M,N,'smallestabs');
w = diag(D).^0.5;
V = VV; 
M_red = diag(diag(V'*M*V)); 
K_red = diag(diag(V'*K*V));
C_red = diag(diag(V'*C*V));
f_red = @(t) V'*f(t);
y0_red = zeros(2*N,1); 
y0_red(1:end/2) = V'*y0(1:end/2);
y0_red(end/2+1:end) = V'*y0(end/2+1:end);

% solution of the problem
[t,y] = ode45(@(t,y) lin_sys(t,y,M_red,C_red,K_red,f_red), tspan, y0_red);

% Recovery with modes acceleration 
K_inv = K^-1;
mm = diag(M_red); 
temp = diag(1./(mm.*(w.^2)));  

for i=1:length(t)
   % Recovery of all the displacements  
   U_modes(:,i) = V*y(i,1:N)' + ((K_inv)-V*temp*V')*f(i);
%    U_direct(:,i) = V*y(i,1:n)';
end


%% Put the displacement in the in.d field 
for j=1:length(model.en)
    for k=1:6
        if model.en(j).c(k) % remove row and column 6(i-1)+k
            index = 6*(j-1)+k;
            U_modes = [U_modes(1:index-1,:);zeros(1,size(U_modes,2));U_modes(index:end,:)]; 
        end
    end            
end

for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*U_modes; 
    for k = 1:model.b(i).nel+1
       model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),end);
%        if we would like to reconstruct the solution over time
%        model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end

end
