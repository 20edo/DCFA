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
%% Generate the wing model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\

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
v = [50:100:1800]; 
q = 1/2*rho.*v.^2; 

% First iteration
[X_old,e_old] = polyeig(K,0*Ca,M);

% Initialize non linear system variables
A=zeros(size(M,1)+1);
b=zeros(size(M,1)+1,1);

% Following iterations
for i=2:length(v)
    
    for k=1:size(X_old,1)
            A(1:size(M,1),1:size(M,1))=e_old(k)^2*M-q(i)/v(i)*e_old(k)*Ca+K-q(i)*Ka;
            A(1:size(M,1),1)=(2*e_old(k)*M-q(i)/v(i)*Ca)*X_old(:,k);
            A(end,1:size(M,1))=2*X_old(:,k)';
            A(end,end)=0;
            b(1:size(M,1),1)=-(e_old(k)^2*M-q(i)/v(i)*e_old(k)*Ca+K-q(i)*Ka)*X_old(:,k);
            b(end)=1-X_old(:,k)'*X_old(:,k);
%             funz=@(z) A*z-b;
%             z0=[X_old(:,k);e_old(k)];
%             z=fsolve(funz,z0);
            z=A\b;
           X(:,k)=X_old(:,k)+z(1:end-1);
           e(k)=e_old(k)+z(end);
    end

    X_old = X;
    e_old = e;
    eig_(i,:) = e_old;
end

