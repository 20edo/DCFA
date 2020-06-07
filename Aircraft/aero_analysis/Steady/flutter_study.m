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
cd aero_analysis\Steady

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
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');

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
v = [50:10:1800]; 
q = 1/2*rho.*v.^2; 

[X_old,e_old] = polyeig(K,0*Ca,M);
% [~, I] = sort(imag(e_old)); 
% e_old = e_old(I); 
% X_old = X_old(:,I); 
% e_old = e_old(n+1:2*n); 
% X_old = X_old(:,n+1:n*2);


for i=2:length(v)
    [X,e] = polyeig(K-q(i)*Ka,-q(i)/v(i)*Ca,M);
%     [~, I] = sort(imag(e));
%     e = e(I);
%     X = X(:,I);
%     e = e(n+1:2*n);
%     X = X(:,n+1:n*2);
    Q=X_old'*X;
    Q_out=Q-diag(diag(Q));
    [~,I]=max(abs(Q));
    X = X(:,I);
    e=e(I);
    X_old = X;
    e_old = e;
    eig_(i,:) = e_old;
end


%% V-g plot 
g = 2*real(eig_)./abs(imag(eig_)); 
figure 
hold on 
for k = 1:n
    plot(v,g(:,k));
end





