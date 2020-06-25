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
%% Performance indicator
% 1 -> Controller based only on modal velocities
% 2 -> Controller based on modal velocities and engines accelerations
% 3 -> Controller based on modal velocities and engines velocities
% 4 -> Controller based on allieviation of the loads at the root of the
%      wing and at the root of the engines' support
controller=1;

switch controller
    
    case 1
        C_z = [zeros(N) eye(N)];
        C_z = C_z/N;
        D_zu = zeros(N,1);
        W_zz = eye(N);
        weight = 0.5;
        weight1 = 0.01;
        weight2 = 0.95;
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
        weight1 = 0.000001;
        weight2 = 0.5;
end

%% Normalizing weight matrixes
% weight=1 ->   Only u counts
% Weight=0 ->   Only z counts

W_zz=(1-weight)*W_zz/norm(W_zz);
W_zz1=(1-weight1)*W_zz/norm(W_zz);
W_zz2=(1-weight2)*W_zz/norm(W_zz);

W_uu = (weight);
W_uu1 = (weight);
W_uu2 = (weight);
%% Setting up the Riccati equation
Q = C_z'*W_zz*C_z;
S = C_z'*W_zz*D_zu;
R = D_zu'*W_zz*D_zu+W_uu;

Q1 = C_z'*W_zz1*C_z;
S1 = C_z'*W_zz1*D_zu;
R1 = D_zu'*W_zz1*D_zu+W_uu1;

Q2 = C_z'*W_zz2*C_z;
S2 = C_z'*W_zz2*D_zu;
R2 = D_zu'*W_zz2*D_zu+W_uu2;

P = are(A-B_u*(R\S'), B_u*(R\B_u'), C_z'*W_zz*C_z-S*(R\S'));
G = R\(B_u'*P+S');

P1 = are(A-B_u*(R1\S1'), B_u*(R1\B_u'), C_z'*W_zz1*C_z-S1*(R1\S1'));
G1 = R1\(B_u'*P1+S1');

P2 = are(A-B_u*(R2\S2'), B_u*(R2\B_u'), C_z'*W_zz2*C_z-S2*(R2\S2'));
G2 = R2\(B_u'*P2+S2');
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

Gain1 = ss([G1 zeros(1,size(C_y,1)-size(G1,2))]);
SYS_controlled1 = feedback(SYS_notcontrolled,Gain1);

Gain2 = ss([G2 zeros(1,size(C_y,1)-size(G2,2))]);
SYS_controlled2 = feedback(SYS_notcontrolled,Gain2);
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
[z1,~,x1] = lsim(SYS_controlled1, u, t);
[z2,~,x2] = lsim(SYS_controlled2, u, t);
[y_nc,~,x_nc] = lsim(SYS_notcontrolled, u, t);
z = z';
z1 = z1';
z2 = z2';
x = x';
y_nc = y_nc';
x_nc = x_nc';
%% I'm recovering the physical coordinates from the modal ones
for i = 1:length(t)
    z_sol(:,i) = -V*z(1:N,i); 
    z_sol_v(:,i) = -V*z(N+1:2*N,i);
    z_sol_acc(:,i) = -V*z(2*N+1:3*N,i);
    
    z_sol1(:,i) = -V*z1(1:N,i); 
    z_sol_v1(:,i) = -V*z1(N+1:2*N,i);
    z_sol_acc1(:,i) = -V*z1(2*N+1:3*N,i);
    
    z_sol2(:,i) = -V*z2(1:N,i); 
    z_sol_v2(:,i) = -V*z2(N+1:2*N,i);
    z_sol_acc2(:,i) = -V*z2(2*N+1:3*N,i);
    
    y_sol_nc(:,i) = -V*y_nc(1:N,i);
    y_sol_nc_v(:,i) = -V*y_nc(N+1:2*N,i);
    y_sol_nc_acc(:,i) = -V*y_nc(2*N+1:3*N,i);
end
z_sol_internalforces= z(3*N+1:end,:);
y_sol_nc_internalforces = y_nc(3*N+1:end,:);

z_sol_internalforces1= z1(3*N+1:end,:);

z_sol_internalforces2= z2(3*N+1:end,:);
%% plot the vertical acceleration of the tip of the wing
% if we want to see the vertical acceleration of the tip of the wing
% (if the model in input is the wing, if not it doesn't have any sense)

if controller == 1
    name = 'velocity';
elseif controller == 4
    name = 'loads';
end
if input == 1 
    ingresso = 'impulse';
elseif input == 2
    ingresso = 'step';
elseif input == 3
    ingresso = 'rect';
end

% if controller == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    plot(t,z_sol(end-3,:))
    hold on
    grid on
    plot(t,z_sol1(end-3,:))
    plot(t,z_sol2(end-3,:))
    plot(t,y_sol_nc(end-3,:),'k')
    title('\quad Vertical displacement of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t \quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$q \quad[m]$','Interpreter','latex','FontSize',14)
    
    subplot(2,2,2)
    plot(t,z_sol_v(end-3,:))
    hold on
    grid on
    plot(t,z_sol_v1(end-3,:))
    plot(t,z_sol_v2(end-3,:))
    plot(t,y_sol_nc_v(end-3,:),'k')
    title('\quad Vertical velocity of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)
    
    subplot(2,2,3)
    plot(t,z_sol_acc(end-3,:))
    hold on
    grid on
    plot(t,z_sol_acc1(end-3,:))
    plot(t,z_sol_acc2(end-3,:))
    plot(t,y_sol_nc_acc(end-3,:),'k')
    title('\quad Vertical acceleration of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,2,4)
    plot(t,-G*x)
    hold on
    grid on
    plot(t,-G1*x)
    plot(t,-G2*x)
    plot(t,u,'k')
    title('Aileron deflection','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\beta\quad[rad]$','Interpreter','latex','FontSize',14)
% elseif controller == 4
%     figure('units','normalized','outerposition',[0 0 1 1])
%     plot(t,-G*x+u)
%     hold on
%     grid on
%     plot(t,-G1*x+u)
%     plot(t,-G2*x+u)
%     plot(t,u,'k')
%     title('Aileron deflection','Interpreter','latex','FontSize',14)
%     xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
%     ylabel('$\beta\quad[rad]$','Interpreter','latex','FontSize',14)
% end

fname = [name '_ideal_tip_' ingresso];
saveas(figure(1),fname,'epsc')

%% Plot the engines' displacements
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
plot(t,z_sol(3,:))
hold on
grid on
plot(t,z_sol1(3,:))
plot(t,z_sol2(3,:))
plot(t,y_sol_nc(3,:),'k')
title('\quad\quad\quad Vertical displacement of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$q$ \quad[m]','Interpreter','latex','FontSize',14)

subplot(2,3,2)
plot(t,z_sol_v(3,:))
hold on
grid on
plot(t,z_sol_v1(3,:))
plot(t,z_sol_v2(3,:))
plot(t,y_sol_nc_v(3,:),'k')
title('\quad\quad Vertical velocity of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)

subplot(2,3,3)
plot(t,z_sol_acc(3,:))
hold on
grid on
plot(t,z_sol_acc1(3,:))
plot(t,z_sol_acc2(3,:))
plot(t,y_sol_nc_acc(3,:),'k')
title('\quad\quad Vertical acceleration of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)


subplot(2,3,4)
plot(t,z_sol(6+3,:))
hold on
grid on
plot(t,z_sol1(6+3,:))
plot(t,z_sol2(6+3,:))
plot(t,y_sol_nc(6+3,:),'k')
title('\quad\quad\quad Vertical displacement of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$q$ \quad[m]','Interpreter','latex','FontSize',14)

subplot(2,3,5)
plot(t,z_sol_v(6+3,:))
hold on
grid on
plot(t,z_sol_v1(6+3,:))
plot(t,z_sol_v2(6+3,:))
plot(t,y_sol_nc_v(6+3,:),'k')
title('\quad\quad Vertical velocity of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)

subplot(2,3,6)
plot(t,z_sol_acc(6+3,:))
hold on
grid on
plot(t,z_sol_acc1(6+3,:))
plot(t,z_sol_acc2(6+3,:))
plot(t,y_sol_nc_acc(6+3,:),'k')
title('\quad\quad Vertical acceleration of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)

fname = [name '_ideal_engine_' ingresso];
saveas(figure(2),fname,'epsc')
%% Plot the stresses for in the CORRENTI for each section
if controller == 4
    % plot of the stress of the 6 CORRENTI in the root section
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(1,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(1,:))
    plot(t,z_sol_internalforces2(1,:))
    plot(t,y_sol_nc_internalforces(1,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(2,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(2,:))
    plot(t,z_sol_internalforces2(2,:))
    plot(t,y_sol_nc_internalforces(2,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(3,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(3,:))
    plot(t,z_sol_internalforces2(3,:))
    plot(t,y_sol_nc_internalforces(3,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(4,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(4,:))
    plot(t,z_sol_internalforces2(4,:))
    plot(t,y_sol_nc_internalforces(4,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(5,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(5,:))
    plot(t,z_sol_internalforces2(5,:))
    plot(t,y_sol_nc_internalforces(5,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(6,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(6,:))
    plot(t,z_sol_internalforces2(6,:))
    plot(t,y_sol_nc_internalforces(6,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_ideal_s1_' ingresso];
    saveas(figure(3),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section before the first
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(7,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(7,:))
    plot(t,z_sol_internalforces2(7,:))
    plot(t,y_sol_nc_internalforces(7,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(8,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(8,:))
    plot(t,z_sol_internalforces2(8,:))
    plot(t,y_sol_nc_internalforces(8,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(9,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(9,:))
    plot(t,z_sol_internalforces2(9,:))
    plot(t,y_sol_nc_internalforces(9,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(10,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(10,:))
    plot(t,z_sol_internalforces2(10,:))
    plot(t,y_sol_nc_internalforces(10,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(11,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(11,:))
    plot(t,z_sol_internalforces2(11,:))
    plot(t,y_sol_nc_internalforces(11,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(12,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(12,:))
    plot(t,z_sol_internalforces2(12,:))
    plot(t,y_sol_nc_internalforces(12,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_ideal_s2_' ingresso];
    saveas(figure(4),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section after the first
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(13,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(13,:))
    plot(t,z_sol_internalforces2(13,:))
    plot(t,y_sol_nc_internalforces(13,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(14,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(14,:))
    plot(t,z_sol_internalforces2(14,:))
    plot(t,y_sol_nc_internalforces(14,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(15,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(15,:))
    plot(t,z_sol_internalforces2(15,:))
    plot(t,y_sol_nc_internalforces(15,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(16,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(16,:))
    plot(t,z_sol_internalforces2(16,:))
    plot(t,y_sol_nc_internalforces(16,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(17,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(17,:))
    plot(t,z_sol_internalforces2(17,:))
    plot(t,y_sol_nc_internalforces(17,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(18,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(18,:))
    plot(t,z_sol_internalforces2(18,:))
    plot(t,y_sol_nc_internalforces(18,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_ideal_s3_' ingresso];
    saveas(figure(5),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section before the second
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(19,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(19,:))
    plot(t,z_sol_internalforces2(19,:))
    plot(t,y_sol_nc_internalforces(19,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(20,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(20,:))
    plot(t,z_sol_internalforces2(20,:))
    plot(t,y_sol_nc_internalforces(20,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(21,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(21,:))
    plot(t,z_sol_internalforces2(21,:))
    plot(t,y_sol_nc_internalforces(21,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(22,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(22,:))
    plot(t,z_sol_internalforces2(22,:))
    plot(t,y_sol_nc_internalforces(22,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(23,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(23,:))
    plot(t,z_sol_internalforces2(23,:))
    plot(t,y_sol_nc_internalforces(23,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(24,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(24,:))
    plot(t,z_sol_internalforces2(24,:))
    plot(t,y_sol_nc_internalforces(24,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_ideal_s4_' ingresso];
    saveas(figure(6),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section after the second
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(25,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(25,:))
    plot(t,z_sol_internalforces2(25,:))
    plot(t,y_sol_nc_internalforces(25,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(26,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(26,:))
    plot(t,z_sol_internalforces2(26,:))
    plot(t,y_sol_nc_internalforces(26,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(27,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(27,:))
    plot(t,z_sol_internalforces2(27,:))
    plot(t,y_sol_nc_internalforces(27,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(28,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(28,:))
    plot(t,z_sol_internalforces2(28,:))
    plot(t,y_sol_nc_internalforces(28,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(29,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(29,:))
    plot(t,z_sol_internalforces2(29,:))
    plot(t,y_sol_nc_internalforces(29,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(30,:))
    hold on
    grid on
    plot(t,z_sol_internalforces1(30,:))
    plot(t,z_sol_internalforces2(30,:))
    plot(t,y_sol_nc_internalforces(30,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_ideal_s5_' ingresso];
    saveas(figure(7),fname,'epsc')
end