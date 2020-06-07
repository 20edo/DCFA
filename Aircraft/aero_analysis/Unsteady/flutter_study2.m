% This is a script to study the flutter problem at 10.000 m for the clamped
% wing

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%
%
%
clear all , close all, clc
cd ..
cd ..
%% Generate the wing model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\Unsteady

% switch off the aerodynamic properties of the engine support
for i=16:19
    aircraft.b(i).ssh = false;
end

%% Build the swept wing model
wing=m_init();
wing.en=[en_ground(aircraft.en(7).x) ...
    aircraft.en(17) aircraft.en(18)];
wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9) aircraft.b(16) aircraft.b(17)];
% wing.en=[en_ground(aircraft.en(7).x)];
% wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9)];
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end

chord=7.72;

wing = m_compute_matrices(wing); 
%% Reduction of the model using n eigenvectors
n = 15;
[V,D] = eigs(wing.K,wing.M,n,'smallestabs');
V_red = V;
M = V'*wing.M*V;
K = V'*wing.K*V;
Ka = V'*wing.Ka*V;
Ca = V'*wing.Ca*V;

%% Altitude fixed to 10.000 m
[T,a,P,rho] = atmosisa(10000);

% the problem is in the form
% M*q_dotdot - q/Vinf*C*q_dot + (K - q*Ka)*q = 0
% the solution of the problem is given by polyeig(K,C,M)

%% Tracking of eigenvalues trough eigenvectors
v = [0:20:400];
q = 1/2*rho.*v.^2;

alpha = 0.01;
gamma = 0.01; 
Cs = alpha*M + gamma*K; 
% First iteration
[X_old,e_old] = polyeig(K-q(1)*Ka,Cs,M);

% Initialize non linear system variables
A=zeros(size(M,1)+1);
b=zeros(size(M,1)+1,1);

% For each airspeed
for i=2:length(v)
    % For each eigenvalue
    for j=1:length(e_old)
        k=imag(e_old(j))*chord/v(i);
        % Assemble Ham matrices
        wing=m_add_unsteady_loads(wing,[v(i) 0 0]',k);
        w=linspace(-1,1,10);
        for k=1:length(w)
        W_vect(:,:,k)=funz(wing,[v(i) 0 0]',w(k));
        end
        Ham_k=gradient(W_vect);
    end
            
    
end

%% V-g plot 
g = 2*real(eig_)./abs(imag(eig_));
figure 
hold on 
for k = 1:n
    plot(v,g(:,k));
end
ylim([-50,50])
ylabel('g')

figure 
hold on 
for k = 1:n
    plot(v,real(eig_(:,k)));
end
ylim([-50,50])
ylabel('real(eig)')

