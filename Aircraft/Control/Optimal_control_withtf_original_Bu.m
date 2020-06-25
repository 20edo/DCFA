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
max_freq = freq(end)/2/pi;
% deltat = 1/(max_freq*8); %Nyquist
deltat = 0.000244141;
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
%% Actuator transfer function
w_cut1=2*pi*5;
w_cut2=2*pi*15;
denominator=conv([1 w_cut1],[1 w_cut2]);
mechanical_actuator=tf(w_cut1*w_cut2, denominator);

figure
bode(mechanical_actuator)
title('Mechanical actuator')
grid on
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
controller=4;

switch controller
    
    case 1
        C_z = [zeros(N) eye(N)];
        D_zu = zeros(N,1);
        risp = abs(freqresp(mechanical_actuator,w));
        risp = risp(1,:);
        W_zz = eye(N).*diag(risp);
        weight = 0.99;
    case 2
        % Define relative importance of modes velocity engines' acceleration
        % 1 -> Only modes velocity count
        % 0 -> Only engines' acceleration count
        rho=1;
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
        %         % Define input matrix in physical coordinates      
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
        % matrix of the physical internal stresses in the considered
        % sections
        C_z=elements_correnti*wing.navier(:,7:end)*V;
        C_z=[C_z zeros(size(C_z))];
        C_z = C_z/sqrt(norm(full(C_z*C_z')));
        D_zu = zeros(size(C_z,1),1);
%         W_zz=eye(size(C_z,1))/(size(C_z,1));
        W_zz = eye(size(C_z,1));
        W_zz(1:3:end,1:3:end) = W_zz(1:3:end,1:3:end)/10;
        weight = 0.01;
end

%% Normalizing weight matrixes
% weight=1 ->   Only u counts
% Weight=0 ->   Only z counts
W_zz=(1-weight)*W_zz/norm(W_zz);

W_uu = (weight);
%% Actuated system
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
SYS_notcontrolled = series(mechanical_actuator,SYS_notcontrolled);
% Add beta and betadot to to the output
C=[SYS_notcontrolled.C;
    zeros(2,size(SYS_notcontrolled.C,2)-2) eye(2)];
D=[SYS_notcontrolled.D; 0; 0];
C_z = [C_z zeros(size(C_z,1),2)];

SYS_notcontrolled = ss(SYS_notcontrolled.A,SYS_notcontrolled.B, C, D);
%% Setting up the Riccati equation
A=SYS_notcontrolled.A;
B_u=SYS_notcontrolled.B;

Q = C_z'*W_zz*C_z;
S = C_z'*W_zz*D_zu;
R = D_zu'*W_zz*D_zu+W_uu;

P = are(A-B_u*inv(R)*S', B_u*inv(R)*B_u', C_z'*W_zz*C_z-S*inv(R)*S');
G = inv(R)*(B_u'*P+S');

SYS_controlled = ss(A-B_u*G,B_u,C,D);
%% Margine di guadagno e di fase di beta
A_beta = SYS_controlled.A;
B_beta = SYS_controlled.B;
SYS = ss(A_beta,B_beta,[zeros(1,20) 1 0], 0);
figure
bode(SYS)
title('Bode of controlled system')
[gm,pm,wcg,wcp] = margin(SYS)
Gain_beta = SYS_controlled.A(N+1:2*N,end)./(M_red\(q*Fb)); %
Gain_beta = abs(Gain_beta(end));
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
%% Update B matrix (the disturbance does not come from the aileron and do not need an actuator
SYS_controlled.B=[zeros(N,1)
    M_red\(q*Fb)
    0
    0];
SYS_notcontrolled.B=[zeros(N,1)
    M_red\(q*Fb)
    0
    0];
%% Plot of the output
q0 = zeros(2*N+2,1);                 
[z,~,x] = lsim(SYS_controlled, u, t); 
[y_nc,~,x_nc] = lsim(SYS_notcontrolled, u, t);
z = z';
x = x';
y_nc = y_nc';
x_nc = x_nc';

%% I'm recovering the physical coordinates from the modal ones
for i = 1:length(t)
    z_sol(:,i) = -V*z(1:N,i); 
    z_sol_v(:,i) = -V*z(N+1:2*N,i);
    z_sol_acc(:,i) = -V*z(2*N+1:3*N,i);
    y_sol_nc(:,i) = -V*y_nc(1:N,i);
    y_sol_nc_v(:,i) = -V*y_nc(N+1:2*N,i);
    y_sol_nc_acc(:,i) = -V*y_nc(2*N+1:3*N,i);
end
z_sol_internalforces= z(3*N+1:end,:);
y_sol_nc_internalforces = y_nc(3*N+1:end,:);
%% beta
[Beta] = lsim(mechanical_actuator,-G*x,t);
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
    plot(t,y_sol_nc(end-3,:),'k')
    title('\quad Vertical displacement of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t \quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$q \quad[m]$','Interpreter','latex','FontSize',14)
    
    subplot(2,2,2)
    plot(t,z_sol_v(end-3,:))
    hold on
    grid on
    plot(t,y_sol_nc_v(end-3,:),'k')
    title('\quad Vertical velocity of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)
    
    subplot(2,2,3)
    plot(t,z_sol_acc(end-3,:))
    hold on
    grid on
    plot(t,y_sol_nc_acc(end-3,:),'k')
    title('\quad Vertical acceleration of the tip','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,2,4)
%    plot(t,x(end,:)*Gain_beta)
    plot(t,Beta)
    hold on
    grid on
    plot(t,u,'k')
    title('Aileron deflection','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\beta\quad[rad]$','Interpreter','latex','FontSize',14)
% elseif controller == 4
%     figure('units','normalized','outerposition',[0 0 1 1])
%     plot(t,Beta)
%     hold on
%     grid on
%     plot(t,u,'k')
%     title('Aileron deflection','Interpreter','latex','FontSize',14)
%     xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
%     ylabel('$\beta\quad[rad]$','Interpreter','latex','FontSize',14)
% end

fname = [name '_tf_extinput_tip_' ingresso];
saveas(figure(3),fname,'epsc')

%% Plot the engines' displacements
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
plot(t,z_sol(3,:))
hold on
grid on
plot(t,y_sol_nc(3,:),'k')
title('\quad\quad\quad\quad Vertical displacement of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$q$ \quad[m]','Interpreter','latex','FontSize',14)

subplot(2,3,2)
plot(t,z_sol_v(3,:))
hold on
grid on
plot(t,y_sol_nc_v(3,:),'k')
title('\quad Vertical velocity of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)

subplot(2,3,3)
plot(t,z_sol_acc(3,:))
hold on
grid on
plot(t,y_sol_nc_acc(3,:),'k')
title('\quad Vertical acceleration of the inner engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)


subplot(2,3,4)
plot(t,z_sol(6+3,:))
hold on
grid on
plot(t,y_sol_nc(6+3,:),'k')
title('\quad Vertical displacement of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$q$ \quad[m]','Interpreter','latex','FontSize',14)

subplot(2,3,5)
plot(t,z_sol_v(6+3,:))
hold on
grid on
plot(t,y_sol_nc_v(6+3,:),'k')
title('\quad Vertical velocity of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\dot{q} \quad[m/s]$','Interpreter','latex','FontSize',14)

subplot(2,3,6)
plot(t,z_sol_acc(6+3,:))
hold on
grid on
plot(t,y_sol_nc_acc(6+3,:),'k')
title('\quad Vertical acceleration of the outer engine','Interpreter','latex','FontSize',14)
xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
ylabel('$\ddot{q} \quad[m/s^{2}]$','Interpreter','latex','FontSize',14)

fname = [name '_tf_extinput_engine_' ingresso];
saveas(figure(4),fname,'epsc')
%% Plot the stresses for in the CORRENTI for each section
if controller == 4
    % plot of the stress of the 6 CORRENTI in the root section
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(1,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(1,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(2,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(2,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(3,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(3,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(4,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(4,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(5,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(5,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(6,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(6,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_tf_extinput_s1_' ingresso];
    saveas(figure(5),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section before the first
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(7,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(7,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(8,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(8,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(9,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(9,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(10,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(10,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(11,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(11,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(12,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(12,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_tf_extinput_s2_' ingresso];
    saveas(figure(6),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section after the first
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(13,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(13,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(14,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(14,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(15,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(15,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(16,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(16,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(17,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(17,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(18,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(18,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_tf_extinput_s3_' ingresso];
    saveas(figure(7),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section before the second
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(19,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(19,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(20,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(20,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(21,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(21,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(22,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(22,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(23,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(23,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(24,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(24,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_tf_extinput_s4_' ingresso];
    saveas(figure(8),fname,'epsc')
    
    % plot of the stress in the 6 CORRENTI in the section after the second
    % engine
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
    plot(t,z_sol_internalforces(25,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(25,:),'k')
    title('Stress 1','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,2)
    plot(t,z_sol_internalforces(26,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(26,:),'k')
    title('Stress 2','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,3)
    plot(t,z_sol_internalforces(27,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(27,:),'k')
    title('Stress 3','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    
    subplot(2,3,4)
    plot(t,z_sol_internalforces(28,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(28,:),'k')
    title('Stress 4','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,5)
    plot(t,z_sol_internalforces(29,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(29,:),'k')
    title('Stress 5','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    subplot(2,3,6)
    plot(t,z_sol_internalforces(30,:))
    hold on
    grid on
    plot(t,y_sol_nc_internalforces(30,:),'k')
    title('Stress 6','Interpreter','latex','FontSize',14)
    xlabel('$t\quad[s]$','Interpreter','latex','FontSize',14)
    ylabel('$\sigma \quad[Pa]$','Interpreter','latex','FontSize',14)
    
    fname = [name '_tf_extinput_s5_' ingresso];
    saveas(figure(9),fname,'epsc')
end