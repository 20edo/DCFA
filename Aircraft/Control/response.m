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
clear all, close all, clc
cd ..
%% Generate the wing model
cd generate_model
generate_model
% Move to the control folder
cd Control\

% switch off the aerodynamic properties of the engine support
for i=16:19 
    aircraft.b(i).ssh = false; 
end

%% Build the swept wing model
wing=m_init();
wing.en=[en_ground(aircraft.en(7).x) ...
    aircraft.en(19) aircraft.en(20)];
wing_list=[aircraft.b(10) aircraft.b(11) aircraft.b(12) aircraft.b(18) aircraft.b(19)];
% wing.en=[en_ground(aircraft.en(7).x)];
% wing_list=[aircraft.b(10) aircraft.b(11) aircraft.b(12)];
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end
%% Calculations for the plotting VTAS and q
[T, a, P, rho] = atmosisa(10000);
Ma = 0.7;
v = Ma*a;
q = 0.5*rho*v^2;
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');
%% Compute properties

wing=m_compute_matrices(wing);

%% Clamp the aircraft (to eliminate rigid body motion)
K = wing.K;
C = wing.C;
M = wing.M;

Ka = wing.Ka;
Ca = wing.Ca;
fb = wing.fb;
fa = wing.fa;

%% Reduced order model
N = 10; % dimension of our reduced model
alpha = 0.3; % shift value
[V, lambda] = eigs(K+alpha*M, M, N,'smallestabs');
% Reduced matrixes
K_red = V'*(K-q*Ka)*V;
C_red = V'*(C-q/v*Ca)*V;
M_red = V'*M*V;

%% Taking into account the max frequency in order to choose a correct 
%  sampling frequency in time
w = sqrt(diag(lambda)-alpha); %prendi doppio freq più alta 
max_freq = w(end);
deltat = 1/(max_freq*2); %Nyquist

%% External input
% I'm taking into account the force given by the aileron deflection
Fb = transpose(V)*fb;

%% Matrix assembly in state space form
A_s = [zeros(N) eye(N);
    -inv(M_red)*K_red -inv(M_red)*C_red];
B_s = [zeros(N,1);
    inv(M_red)*q*Fb];
C_s = [diag(ones(N,1)) zeros(N);
    zeros(N) zeros(N)];
D_s = 0;

%% State space system
SYS = ss(A_s,B_s,C_s,D_s);
% Now we want the transfer function
SYS_tf = tf(SYS);

%% Number of problem to be solved
problem=3;
%% #1
% sinusoidal input
%% #2
% impulse input
%% #3
% step input
%% #4 
% aileron deflection beta for prescribed asymptotic roll rate 
%% #5 
% initial roll acc p_dot for prescribed roll rate p
%% 6 
% roll rate p required for prescribed roll acceleration 

%% Define the intput & output
% t = [0:deltat:5];
% al posto di 10 fai 2pi/w(1)
t = [0:0.1:10];

switch problem
    case 1
        u=sin(t)'; 
        [y]=lsim(SYS,u,t);
        
    case 2
        [y] = impulse(SYS_tf,t(end));
        
    case 3
        [y] = step(SYS_tf,t);
    case 4 
        %% free the pitch dof
        % wing.en(1).c(5) = false;
    case 5 
        
    case 6 
           
end

for i = 1:length(t)
    y_sol(:,i) = V*y(i,1:N)'; % I'm recovering the physic coordinates from
                              % the modal ones
end
y_sol = [y_sol];
y_sol_static = [(K-q*Ka)\(q*fa)];
%% plot the vertical displacement of the tip of the wing
% if we want to see the vertical displacement of the tip of the wing 
% (if the model in input is the wing, if not it doesn't have any sense)
% u_tot = u_sol+u_sol_static;
% plot(t,u_tot(end-3,:))
% xlabel('t(s)')
% ylabel('w(m)')

%% plot of the model with only static aero loads
if 0
    for i=4:5
        wing.b(i).ssh = true;
    end
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;         % we have only one eig
    m_plot_eigenshape(wing,options,y_sol_static*10);
    m_plot_eigenshape_easy(wing,y_sol_static*10)
end
%% video of the easy-model in response to the input
if 0
    FR=m_plot_video_easy(wing,t,y_sol*10+y_sol_static*10); %y_sol*100+y_sol_static*10
end




