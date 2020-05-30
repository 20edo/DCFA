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
function model = m_solution_dynamic_free(model,tspan,initialcond,f,N, alpha)

 if (~exist('alpha', 'var'))
     % "dof" parameter does not exist
      alpha = 0.2;
 end
cd .. 
[w, U, k] = ROM_solver(N, model.M, model.K, alpha);
cd Model
M = model.M; 
K = model.K; 
C = model.C; 
y0 = initialcond; 

%% Way relatively fast to compute how many rigid modes we have 
n_node = size(model.en,2);
dof_node = zeros(1,n_node);
for i=1:size(model.b,2)
    a = model.b(i).on; 
    b = model.b(i).en; 
    dof_node(1,a) = dof_node(1,a) + 1;
    dof_node(1,b) = dof_node(1,b) + 1;
end
% compute the number of degrees of freedom of the entire structure
dof  = sum([model.b(:).nel]+1)-sum((dof_node(1,:)-1)); 
dof = 6*dof; 
number_rigid_modes = 6 - (dof - size(model.M,1)); 
N_def = N - number_rigid_modes; % dimention of the deformative problem

%% 
% divide the eigenvalues in rigid and deformative
U_r = U(:,1:number_rigid_modes); 
U_d = U(:,number_rigid_modes+1:end); 

% prepare Projector matrix and K-tilda
K_tilda = K + sum(diag(K))/length(diag(K))*U_r*eye(number_rigid_modes)*U_r'; 
P = eye(dof) - M*U_r*(U_r'*M*U_r)^-1*U_r'; 

% solve the deformative problem 
M_def = diag(diag(U_d'*M*U_d)); 
K_def = diag(diag(U_d'*K*U_d));
C_def = diag(diag(U_d'*C*U_d));
f_def = @(t) U_d'*f(t); 
y0_def = zeros(2*N_def,1); 
y0_def(1:end/2) = U_d'*y0(1:end/2);
y0_def(end/2+1:end) = U_d'*y0(end/2+1:end);

[t,y] = ode45(@(t,y) lin_sys(t,y,M_def,C_def,K_def,f_def), tspan, y0_def);

temp = K_tilda^-1;
A = P'*(temp)*P;
B = P'*(temp)*M*U_d;
for i=1:length(t)
   q_d_dot_dot(:,i) = 1./diag(M_def).*(-diag(C_def).*y(i,N_def+1:end)'-diag(K_def).*y(i,1:N_def)'+f_def(i));
   u(:,i) = A*f(i) - B*q_d_dot_dot(:,i); 
end

for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*u; 
    for k = 1:model.b(i).nel+1
%        model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),end);
%        if we would like to reconstruct the solution over time
       model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end
