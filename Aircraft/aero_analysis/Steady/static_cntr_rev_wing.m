% Static aero analysis of the clamped wing
% - Control reversal
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
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end

%% Generate the straight wing model
wing_straight = wing;
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');
wing_straight = m_add_aero_loads_straight(wing_straight,[1,0,0]');


% wing.b(1).fb = zeros(size(wing.b(1).fb));
% wing.b(1).Lq = zeros(size(wing.b(1).Lq));
% wing.b(1).Lb = zeros(size(wing.b(1).Lb));
% wing.b(2).fb = zeros(size(wing.b(2).fb));
% wing.b(2).Lq = zeros(size(wing.b(2).Lq));
% wing.b(2).Lb = zeros(size(wing.b(2).Lb));
%
% wing = m_compute_matrices(wing);



%% Assebly stiffness the K matrices for cntr_rev
K_cr = [wing.K, zeros(size(wing.K,1),1); zeros(1,size(wing.K,2)),0];
Ka_cr = [wing.Ka, wing.fb; wing.Lq, wing.Lb];

%% Direct problem
q = 24500;
A = [wing.K-q*wing.Ka, zeros(size(wing.K,1),1); -q*wing.Lq, 1];
b = q*[wing.fb; wing.Lb]*deg2rad(20);
sol = A\b;
DL = sol(end);


%% Find the control reversal (cr) dynamic pressure
[V_cr,D_cr]= eig(full(K_cr),full(Ka_cr));
q_cr = diag(D_cr);
q_cr = q_cr.*(abs(imag(q_cr))<10^-3);
[q_cr,I] = sort(real(q_cr));
V_cr = V_cr(:,I);                         % sort the eigenshapes
V_cr(:,q_cr<1)=[];                    % select the eigenshapes with positive eig
q_cr(q_cr<1)=[];
q_cr = q_cr(1);                   % select the minimum positive q_inf
V_cr = V_cr(:,1);                     % select its eigenshape


%% Calculations for the plotting VTAS and MACH when altitude changes - swept wing
[T, a, P, rho] = atmosisa(0:100:11000);v = sqrt(q_cr*2./rho);
M = v./a;


figure

set(gcf, 'Position',  [40, 40, 700, 500])
subplot(1,2,1)
plot(0:100:11000,v,'LineWidth',2)
grid on
xlabel('Altitude [m]')
ylabel('Divergence Velocity [m/s]')
title('Divergence Velocity Swept Wing')
subplot(1,2,2)
plot(0:100:11000,M,'LineWidth',2)
grid on
xlabel('Altitude [m]')
ylabel('Divergence Mach [-]')
title('Divergence Mach Swept Wing')

if 1
    % switch on the aero properties for the plot
    for i=4:5
        wing.b(i).ssh = true;
    end
    
    phi = rad2deg(0);
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    m_plot_eigenshape(wing,options,((V_cr(1:end-1)))*3);
end
