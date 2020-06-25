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
% k = linspace(0.1,1e-12,500); 
temp =  1-cos(linspace(0,pi/2,100000000)); 
k = temp(1:150); 
k = flip(k); 
k = k(1:end-1); 
error = zeros(1,length(k)); 
wing=m_add_unsteady_loads(wing,[v(i) 0 0]',0);
Ham_zero = wing.Ham; 
wing=m_add_unsteady_loads(wing,[v(i) 0 0]',k(1));
derivative_old = imag(wing.Ham)/k(1); 
% For each airspeed
for i=2:length(k)
    tic
    % Assemble Ham matrices
    wing=m_add_unsteady_loads(wing,[1 0 0]',k(i));
    derivative = imag(wing.Ham)/k(i);
    error(i) = norm(full(derivative - derivative_old))/norm(full(derivative_old)); 
    derivative_old = derivative;
    disp(i)
    toc
end

plot(k,error,'-o'); 


