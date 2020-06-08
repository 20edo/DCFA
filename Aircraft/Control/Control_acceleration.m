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
alpha = 0;
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

%% Matrix assembly
A = [zeros(N) eye(N);
    -M_red\K_red -M_red\C_red];

B_u = [zeros(N,1);
    M_red\(q*Fb)];

B_d = [zeros(N,1);
    M_red\ones(N,1)];

C_y = [eye(N) zeros(N)];

C_yv = [zeros(N) eye(N)];

D_yu = 0;

D_yd = 0;

%% Performance indicator
% % C_z = [-M_red\K_red -M_red\C_red]; 
% % D_zu = [M_red\(q*Fb)];

C_z = [zeros(N) eye(N)];
D_zu = zeros(N,1);

% C_z = [zeros(N) eye(N)]*A;
% C_z = A(N+1:end,:);
% D_zu = B_u(N+1:end);
%% Weight matrixes
weight = 0.1;

W_zz = (1-weight) * sqrt(lambda)/(sum(sum(sqrt(lambda))));
% W_zz = (1-weight) * (lambda)/(sum(sum(lambda)));

W_uu = (weight);
%% Setting up the Riccati equation
Q = C_z'*W_zz*C_z;
S = C_z'*W_zz*D_zu;
R = D_zu'*W_zz*D_zu+W_uu;

P = are(A-B_u*inv(R)*S', B_u*inv(R)*B_u', C_z'*W_zz*C_z-S*inv(R)*S');
G = inv(R)*(B_u'*P+S');

% Ritar
%% State space model
A_controlled = A-B_u*G;

SYS_controlled = ss(A_controlled, B_u, C_y, D_yu);
SYS_controlled_velocity = ss(A_controlled, B_u, C_yv, D_zu);
SYS_controlled_accelerations = ss(A_controlled, B_u, C_z, D_zu);
SYS_notcontrolled = ss(A, B_u, C_y, D_yu); % Displacements
SYS_notcontrolled_velocity = ss(A, B_u, C_yv, D_zu);
SYS_notcontrolled_accelerations = ss(A, B_u, C_z, D_zu);
%% Plot of the output
t = [0:deltat:10];
beta = deg2rad(5);
% u = [zeros(N,length(t));
%     beta*eye(N,length(t))]; 
% 
% u_func = @(t) beta.*(t<2).*(t>1);

% u_func = @(t) [zeros(N,1);
%     beta*ones(N,1)].*(randn(1));


% u_func = @(t) [zeros(N,1);
%     beta*ones(N,1)].*t;
u_func= @(t) beta.*(t<1);
u=[];
for i =1:length(t)
    u(i)=u_func(t(i));  
    
end

[z,~,x] = lsim(SYS_controlled, u, t);
[z_v] = lsim(SYS_controlled_velocity, u, t);
[z_acc] = lsim(SYS_controlled_accelerations, u, t);
[y_nc,~,x_nc] = lsim(SYS_notcontrolled, u, t);
[y_nc_v] = lsim(SYS_notcontrolled_velocity, u, t);
[y_nc_acc] = lsim(SYS_notcontrolled_accelerations, u, t);

for i = 1:length(t)
    z_sol(:,i) = V*z(i,1:N)'; % I'm recovering the physical coordinates from
                              % the modal ones
    z_sol_v(:,i) = V*z_v(i,1:N)';
    z_sol_acc(:,i) = V*z_acc(i,1:N)';
    y_sol_nc(:,i) = V*y_nc(i,1:N)';
    y_sol_nc_v(:,i) = V*y_nc_v(i,1:N)';
    y_sol_nc_acc(:,i) = V*y_nc_acc(i,1:N)';
end

%% plot the vertical acceleration of the tip of the wing
% if we want to see the vertical acceleration of the tip of the wing 
% (if the model in input is the wing, if not it doesn't have any sense)
figure

subplot(2,2,1)
plot(t,z_sol(end-3,:))
hold on
grid on
plot(t,y_sol_nc(end-3,:))
title('Vertical displacement of the tip')
legend('displacements cotrolled','displacements not controlled')
xlabel('t(s)')
ylabel('$q$','Interpreter','latex')


subplot(2,2,2)
plot(t,z_sol_v(end-3,:))
hold on
grid on
plot(t,y_sol_nc_v(end-3,:))
title('Vertical velocity of the tip')
legend('velocities cotrolled','velocities not controlled')
xlabel('t(s)')
ylabel('$\dot{q}$','Interpreter','latex')

subplot(2,2,3)
plot(t,z_sol_acc(end-3,:))
hold on
grid on
plot(t,y_sol_nc_acc(end-3,:))
title('Vertical acceleration of the tip')
legend('accelerations cotrolled','accelerations not controlled')
xlabel('t(s)')
ylabel('$\ddot{q}$','Interpreter','latex')


subplot(2,2,4)
plot(t,-G*x')
hold on
grid on
plot(t,-G*x'+u)
plot(t,u)
title('Aileron deflection')
legend('Controller output','Controlled','Non controlled')
xlabel('t(s)')
ylabel('$\beta$','Interpreter','latex')


%% Plot the engines' displacements

plot(t,z_sol(6+3,:))
hold on
grid on
plot(t,y_sol_nc(6+3,:))
title('Vertical displacement of the engine')
legend('displacements controlled','displacements not controlled')
xlabel('t(s)')
ylabel('$q$','Interpreter','latex')

plot(t,z_sol(12+3,:))
hold on
grid on
plot(t,y_sol_nc(12+3,:))
title('Vertical displacement of the engine')
legend('displacements controlled','displacements not controlled')
xlabel('t(s)')
ylabel('$q$','Interpreter','latex')


