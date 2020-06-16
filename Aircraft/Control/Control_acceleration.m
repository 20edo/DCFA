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

D_yu = 0;

D_yd = 0;

%% Performance indicator
% % C_z = [-M_red\K_red -M_red\C_red];
% % D_zu = [M_red\(q*Fb)];
% 1 -> Controller based only on modal velocities
% 2 -> Controller based on modal velocities and engines accelerations
% 3 -> Controller based on modal velocities and engines velocities
% 4 -> Controller based on allieviation of the loads at the root of the
%      wing and at the root of the engines' support
controller=4;

switch controller
    
    case 1
        C_z = [zeros(N) eye(N)];
        C_z = C_z/N;
        D_zu = zeros(N,1);
        W_zz = sqrt(lambda)/(sum(sum(sqrt(lambda))));
    case 2
        % Define relative importance of modes velocity engines' acceleration
        % 1 -> Only modes velocity count
        % 0 -> Only engines' acceleration count
        rho=1;
        %
        Minv=inv(M);
        index=1:12;
        % Matrix to value engines' acceleration in physical coordinates
        C_z_phy=[-Minv(index,index)*K(index,index)  zeros(index(end),size(K,1)-index(end))       -Minv(index,index)*C(index,index)    zeros(index(end),size(K,1)-index(end)) ;
            zeros(2*size(K,1)-index(end),2*size(K,1))];
        % Matrix to value the engines' acceleration in reduced coordinates (modal)
        C_z_e=blkdiag(V',V')*C_z_phy*blkdiag(V,V);
        % Eliminate unusefull rows
        C_z_e(N+1:end,:)=[];
        % Matrix to value the velocity of the whole wing
        C_z_v= [zeros(N) eye(N)];
        % Normalize matrices
        C_z_e=C_z_e/norm(C_z_e*C_z_e');
        C_z_v=C_z_v/norm(C_z_v*C_z_v');
        C_z=(1-rho)*C_z_e+rho*C_z_v;
        %         % Define input matrix in physical coordinates      NOT NECESSARY,
        %         THERE IS NO AERO-FORCE ON THE ENGINES
        %         D_zu_phy=[Minv(index,index)*q*fb(index) zeros(index(end),size(K,1)-index(end))];
        % Define the input matrix
        D_zu = zeros(N,1);
        W_zz = sqrt(lambda)/(sum(sum(sqrt(lambda))));
    case 3
        % Define relative importance of modes velocity engines' acceleration
        % 1 -> Only modes velocity count
        % 0 -> Only engines' velocity count
        rho=1;
        %
        Minv=inv(M);
        index=1:12;
        % Matrix to value engines' veclocities in physical coordinates
        C_z_phy=[-zeros(index(end),index(end))  zeros(index(end),size(K,1)-index(end))       eye(index(end))    zeros(index(end),size(K,1)-index(end)) ;
            zeros(2*size(K,1)-index(end),2*size(K,1))];
        % Matrix to value the engines' velocities in reduced coordinates (modal)
        C_z_e=blkdiag(V',V')*C_z_phy*blkdiag(V,V);
        % Eliminate unusefull rows
        C_z_e(N+1:end,:)=[];
        % Matrix to value the velocity of the whole wing
        C_z_v= [zeros(N) eye(N)];
        % Normalize matrices
        C_z_e=C_z_e/norm(C_z_e*C_z_e');
        C_z_v=C_z_v/norm(C_z_v*C_z_v');
        C_z=(1-rho)*C_z_e+rho*C_z_v;
        %         % Define input matrix in physical coordinates      NOT NECESSARY,
        %         THERE IS NO AERO-FORCE ON THE ENGINES
        %         D_zu_phy=[Minv(index,index)*q*fb(index) zeros(index(end),size(K,1)-index(end))];
        % Define the input matrix
        D_zu = zeros(N,1);
        W_zz = sqrt(lambda)/(sum(sum(sqrt(lambda))));
    case 4
        % Select the elements where I want to minimize the stress (root,
        % element before and after the engines' support)
        elements_correnti=[wing.b(1).A(:,[7:12 end-6:end]) wing.b(2).A(:,[7:12 end-6:end]) wing.b(3).A(:,[7:12])];
        elements_correnti=elements_correnti';
        C_z=elements_correnti*wing.navier(:,7:end)*K*V;
        C_z=[C_z zeros(size(C_z))];
%         D_Ku = [eye(N)-M*V*diag(1/M_red)*V'];
        D_zu = zeros(size(C_z,1),1);
        W_zz=eye(size(C_z,1));
end

% C_z = [zeros(N) eye(N)]*A;
% C_z = A(N+1:end,:);
% D_zu = B_u(N+1:end);
%% Normalizing weight matrixes
weight = 0.1;
% weight=1 ->   Only u counts
% Weight=0 ->   Only z counts

W_zz=(1-weight)*W_zz/norm(W_zz);
% W_zz = (1-weight) * (lambda)/(sum(sum(lambda)));

W_uu = (weight);

%% Setting up the Riccati equation
Q = C_z'*W_zz*C_z;
S = C_z'*W_zz*D_zu;
R = D_zu'*W_zz*D_zu+W_uu;

P = are(A-B_u*inv(R)*S', B_u*inv(R)*B_u', C_z'*W_zz*C_z-S*inv(R)*S');
G = inv(R)*(B_u'*P+S');

%% State space model
A_controlled = A-B_u*G;
% Simulate the model to obtain displacements
SYS_controlled = ss(A_controlled, B_u, C_y, D_yu);
SYS_notcontrolled = ss(A, B_u, C_y, D_yu);

%% Recover velocities and accelerations
% Build velocity recover matrix
C_Recover_v = [zeros(N) eye(N)];
D_Recover_v=0;
% Recover velocities
SYS_controlled_velocity = ss(A_controlled, B_u, C_Recover_v, D_Recover_v);
SYS_notcontrolled_velocity = ss(A, B_u, C_Recover_v, D_Recover_v);
% Build acceleration recover matrix
C_Recover_acc = [-M_red\K_red -M_red\C_red];
D_Recover_acc = [M_red\(q*Fb)];
% Recover accelerations
SYS_controlled_accelerations = ss(A_controlled, B_u, C_Recover_acc, D_Recover_acc);
SYS_notcontrolled_accelerations = ss(A, B_u, C_Recover_acc, D_Recover_acc);

%% Recover stresses 
elements_correnti=[wing.b(1).A(:,[7:12 end-6:end]) wing.b(2).A(:,[7:12 end-6:end]) wing.b(3).A(:,[7:12])];
elements_correnti=elements_correnti';
C_Recover_internalforces= elements_correnti*wing.navier(:,7:end)*K*V;
C_Recover_internalforces=[C_Recover_internalforces zeros(size(C_Recover_internalforces))];
D_Recover_internalforces = zeros(size(C_Recover_internalforces,1),1);
SYS_controlled_internalforces = ss(A_controlled, B_u, C_Recover_internalforces, D_Recover_internalforces);
SYS_notcontrolled_internalforces = ss(A, B_u, C_Recover_internalforces, D_Recover_internalforces);

%% Define the input
% Define time axis
t = [0:deltat:1];
beta = deg2rad(2);
% input =
% 1     -> impulse
% 2     -> step
% 3     -> rect
% 4     -> rampa
% 5     -> randn
% 6     -> sin
input=3;
switch input
    case 1
        u_func= @(t) beta.*(t<=0);
    case 2
        u_func= @(t) beta;
    case 3
        u_func= @(t) beta.*(t<=1);
    case 4
        u_func= @(t) beta.*t;
    case 5
        rng(0);
        u_func= @(t) beta.*randn(1);
    case 6
        u_func= @(t) beta.*sin(20*t);
    case 7
        freq=100;
        u_func= @(t) beta.*sinc(freq*t-freq/2);
end

u=[];
for i =1:length(t)
    u(i)=u_func(t(i));
    
end
%% Plot of the output
[z,~,x] = lsim(SYS_controlled, u, t);
[z_v] = lsim(SYS_controlled_velocity, u, t);
[z_acc] = lsim(SYS_controlled_accelerations, u, t);
[z_Recover_internalforces] = lsim(SYS_controlled_Recover_internalforces, u, t);
[y_nc,~,x_nc] = lsim(SYS_notcontrolled, u, t);
[y_nc_v] = lsim(SYS_notcontrolled_velocity, u, t);
[y_nc_acc] = lsim(SYS_notcontrolled_accelerations, u, t);
[y_nc_Recover_internalforces] = lsim(SYS_notcontrolled_Recover_internalforces, u, t);

for i = 1:length(t)
    z_sol(:,i) = V*z(i,1:N)'; % I'm recovering the physical coordinates from
    % the modal ones
    z_sol_v(:,i) = V*z_v(i,1:N)';
    z_sol_acc(:,i) = V*z_acc(i,1:N)';
    z_sol_Recover_internalforces(:,i) = V*z_Recover_internalforces(i,1:N)';
    y_sol_nc(:,i) = V*y_nc(i,1:N)';
    y_sol_nc_v(:,i) = V*y_nc_v(i,1:N)';
    y_sol_nc_acc(:,i) = V*y_nc_acc(i,1:N)';
    y_sol_nc_Recover_internalforces(:,i) = V*y_nc_Recover_internalforces(i,1:N)';
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

figure
subplot(2,3,1)
plot(t,z_sol(3,:))
hold on
grid on
plot(t,y_sol_nc(3,:))
title('Vertical displacement of the inner engine')
legend('displacements controlled','displacements not controlled')
xlabel('t(s)')
ylabel('$q$','Interpreter','latex')

subplot(2,3,2)
plot(t,z_sol_v(3,:))
hold on
grid on
plot(t,y_sol_nc_v(3,:))
title('Vertical velocity of the inner engine')
legend('Velocity controlled','Velocity not controlled')
xlabel('t(s)')
ylabel('$\dot{q}$','Interpreter','latex')

subplot(2,3,3)
plot(t,z_sol_acc(3,:))
hold on
grid on
plot(t,y_sol_nc_acc(3,:))
title('Vertical acceleration of the inner engine')
legend('Acceleration controlled','Acceleration not controlled')
xlabel('t(s)')
ylabel('$\ddot{q}$','Interpreter','latex')


subplot(2,3,4)
plot(t,z_sol(6+3,:))
hold on
grid on
plot(t,y_sol_nc(6+3,:))
title('Vertical displacement of the outer engine')
legend('displacements controlled','displacements not controlled')
xlabel('t(s)')
ylabel('$q$','Interpreter','latex')

subplot(2,3,5)
plot(t,z_sol_v(6+3,:))
hold on
grid on
plot(t,y_sol_nc_v(6+3,:))
title('Vertical velocity of the outer engine')
legend('Velocity controlled','Velocity not controlled')
xlabel('t(s)')
ylabel('$\dot{q}$','Interpreter','latex')

subplot(2,3,6)
plot(t,z_sol_acc(6+3,:))
hold on
grid on
plot(t,y_sol_nc_acc(6+3,:))
title('Vertical acceleration of the outer engine')
legend('Acceleration controlled','Acceleration not controlled')
xlabel('t(s)')
ylabel('$\ddot{q}$','Interpreter','latex')

%% Plot the internal forces at the root

