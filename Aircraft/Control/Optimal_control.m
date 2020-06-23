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
%% Add weight loads
wing = m_add_mass_forces(wing, 1, [0 0 -1]');
%% Compute properties

wing=m_compute_matrices(wing);
% wing.navier = eye(size(wing.navier));
%% Matrix of the model
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
%% Taking into account the max frequency in order to choose a correct
%  sampling frequency in time
w = sqrt(diag(lambda)-alpha); 
freq = w/2/pi;
max_w = w(end);
max_freq = max_w/2/pi;
% deltat = 1/(max_freq*8); %Nyquist
deltat = 0.001;
%% Reduced matrixes
K_red = V'*(K-q*Ka)*V;
C_red = V'*(C-q/v*Ca)*V;
M_red = V'*M*V;

% csi=2e-2;
% C_struct = 2*csi*diag(freq'*M_red);
% C_red = C_red + C_struct;
%% External input
% I'm taking into account the force given by the aileron deflection
Fb = transpose(V)*fb;

%% Matrix assembly of the model
A = [zeros(N) eye(N);
    -M_red\K_red -M_red\C_red];

B_u = [zeros(N,1);
    M_red\(q*Fb)];

% B_d = [zeros(N,1);
%     M_red\ones(N,1)];

%% Performance indicator
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
        weight = 0.5;
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
        C_z_e=C_z_e/sqrt(norm(C_z_e*C_z_e'));
        C_z_v=C_z_v/sqrt(norm(C_z_v*C_z_v'));
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
        C_z_e=C_z_e/sqrt(norm(C_z_e*C_z_e'));
        C_z_v=C_z_v/sqrt(norm(C_z_v*C_z_v'));
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
        elements_correnti=[wing.b(1).A(:,[7:12 end-11:end-6]) wing.b(2).A(:,[7:12 end-11:end-6]) wing.b(3).A(:,[7:12])];
        elements_correnti=elements_correnti';
        C_z=elements_correnti*wing.navier(:,7:end)*V;
        C_z=[C_z zeros(size(C_z))];
        C_z = C_z/sqrt(norm(full(C_z*C_z')));
        D_zu = zeros(size(C_z,1),1);
        W_zz=eye(size(C_z,1))/(size(C_z,1));
        weight = 0.01;
end

%% Normalizing weight matrixes
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
  
%% Recover displacements, velocities, accelerations and physical stresses 
elements_correnti=[wing.b(1).A(:,[7:12 end-11:end-6]) wing.b(2).A(:,[7:12 end-11:end-6]) wing.b(3).A(:,[7:12])];
elements_correnti=elements_correnti';
C_stresses=elements_correnti*wing.navier(:,7:end)*V;

C_y = [eye(N) zeros(N); %displacements
    zeros(N) eye(N); %velocities
    -M_red\K_red -M_red\C_red; %accelerations
    C_stresses zeros(size(C_stresses))]; %internal stresses

D_yu = [zeros(N,1);
    zeros(N,1);
    M_red\(q*Fb);
    zeros(size(C_stresses,1),1)];

SYS_notcontrolled = ss(A, B_u, C_y, D_yu);
% SYS_controlled = ss(A_controlled, B_u, C_y, D_yu);
Gain = ss([G zeros(1,size(C_y,1)-size(G,2))]);
SYS_controlled = feedback(SYS_notcontrolled,Gain);

%% Define the input
% Define time axis
t = [0:deltat:7];
beta = deg2rad(2);
% input =
% 1     -> impulse
% 2     -> step
% 3     -> rect
% 4     -> rampa
% 5     -> randn
% 6     -> sin
% 7     -> sinc
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
        freq=2;
        u_func= @(t) beta.*sinc(freq*t-freq/2);
end

u=[];
for i =1:length(t)
    u(i)=u_func(t(i));
    
end
%% Plot of the output                   
[z,~,x] = lsim(SYS_controlled, u, t);
[y_nc,~,x_nc] = lsim(SYS_notcontrolled, u, t);
z = z';
x = x';
y_nc = y_nc';
x_nc = x_nc';
%% I'm recovering the physical coordinates from the modal ones
for i = 1:length(t)
    z_sol(:,i) = V*z(1:N,i); 
    z_sol_v(:,i) = V*z(N+1:2*N,i);
    z_sol_acc(:,i) = V*z(2*N+1:3*N,i);
    y_sol_nc(:,i) = V*y_nc(1:N,i);
    y_sol_nc_v(:,i) = V*y_nc(N+1:2*N,i);
    y_sol_nc_acc(:,i) = V*y_nc(2*N+1:3*N,i);
end
z_sol_internalforces= z(3*N+1:end,:);
y_sol_nc_internalforces = y_nc(3*N+1:end,:);

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
legend('displacements cotrolled','displacements not controlled','Location','southeast')
xlabel('t[s]')
ylabel('$q$','Interpreter','latex')


subplot(2,2,2)
plot(t,z_sol_v(end-3,:))
hold on
grid on
plot(t,y_sol_nc_v(end-3,:))
title('Vertical velocity of the tip')
legend('velocities cotrolled','velocities not controlled')
xlabel('$t\quad[s]$','Interpreter','latex')
ylabel('$\dot{q}\quad[m]$','Interpreter','latex')

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
plot(t,-G*x)
hold on
grid on
plot(t,-G*x+u)
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

%% Plot the stresses for in the CORRENTI for each section
if controller == 4
    % plot of the stress of the 6 CORRENTI in the root section
    figure
    subplot(2,3,1)
    plot(t,z_sol_internalforces(1,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(1,:))
    title('Stress 1 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(2,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(2,:))
    title('Stress 2 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(3,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(3,:))
    title('Stress 3 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(4,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(4,:))
    title('Stress 4 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(5,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(5,:))
    title('Stress 5 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(6,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(6,:))
    title('Stress 6 at root section')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    % plot of the stress in the 6 CORRENTI in the section before the first
    % engine
    figure
    subplot(2,3,1)
    plot(t,z_sol_internalforces(7,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(7,:))
    title('Stress 1 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(8,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(8,:))
    title('Stress 2 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(9,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(9,:))
    title('Stress 3 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(10,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(10,:))
    title('Stress 4 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(11,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(11,:))
    title('Stress 5 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(12,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(12,:))
    title('Stress 6 at the section before the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    % plot of the stress in the 6 CORRENTI in the section after the first
    % engine
    figure
    subplot(2,3,1)
    plot(t,z_sol_internalforces(13,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(13,:))
    title('Stress 1 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(14,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(14,:))
    title('Stress 2 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(15,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(15,:))
    title('Stress 3 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(16,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(16,:))
    title('Stress 4 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(17,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(17,:))
    title('Stress 5 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(18,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(18,:))
    title('Stress 6 at the section after the first engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    % plot of the stress in the 6 CORRENTI in the section before the second
    % engine
    figure
    subplot(2,3,1)
    plot(t,z_sol_internalforces(19,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(19,:))
    title('Stress 1 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(20,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(20,:))
    title('Stress 2 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(21,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(21,:))
    title('Stress 3 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(22,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(22,:))
    title('Stress 4 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(23,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(23,:))
    title('Stress 5 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(24,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(24,:))
    title('Stress 6 at the section before the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    % plot of the stress in the 6 CORRENTI in the section after the second
    % engine
    figure
    subplot(2,3,1)
    plot(t,z_sol_internalforces(25,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(25,:))
    title('Stress 1 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(26,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(26,:))
    title('Stress 2 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(27,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(27,:))
    title('Stress 3 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(28,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(28,:))
    title('Stress 4 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(29,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(29,:))
    title('Stress 5 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(30,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(30,:))
    title('Stress 6 at the section after the second engine')
    legend('stress controlled','stress not controlled')
    xlabel('t(s)')
    ylabel('$\sigma$','Interpreter','latex')
end