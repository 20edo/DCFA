% Static aero analysis of the T-Tail
% - Control Reversal
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
%% Generate the aircraft model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\Steady


%% Build the wing model
ttail=m_init();
ttail.en=[en_ground(aircraft.en(5).x)];
ttail_list=[aircraft.b(13) aircraft.b(14) aircraft.b(15)];
for i=1:length(ttail_list)
    ttail=m_add_beam(ttail,ttail_list(i));
end

%% Add the aero loads
ttail = m_add_aero_loads(ttail,[1,0,0]');

%% Control reversal for the rudder
% clear fb, Lq, Lb of the stabilizer
rudder = ttail;
rudder.b(2).fb = zeros(size(rudder.b(2).fb));
rudder.b(2).Lq = zeros(size(rudder.b(2).Lq));
rudder.b(2).Lb = zeros(size(rudder.b(2).Lb));
rudder.b(3).fb = zeros(size(rudder.b(3).fb));
rudder.b(3).Lq = zeros(size(rudder.b(3).Lq));
rudder.b(3).Lb = zeros(size(rudder.b(3).Lb));

rudder = m_compute_matrices(rudder);

% Assebly stiffness the K matrices for cntr_rev
K_cr_r = [rudder.K, zeros(size(rudder.K,1),1); zeros(1,size(rudder.K,2)),0];
Ka_cr_r = [rudder.Ka, rudder.fb; rudder.Lq, rudder.Lb];

% Find the control reversal (cr) dynamic pressure
[V_cr_r,D_cr_r]= eig(full(K_cr_r),full(Ka_cr_r));
q_cr_r = diag(D_cr_r);
q_cr_r = q_cr_r.*(abs(imag(q_cr_r))<10^-3);
[q_cr_r,I] = sort(real(q_cr_r));
V_cr_r = V_cr_r(:,I);                         % sort the eigenshapes
V_cr_r(:,q_cr_r<1)=[];                    % select the eigenshapes with positive eig
q_cr_r(q_cr_r<1)=[];
q_cr_r = q_cr_r(1);                   % select the minimum positive q_inf
V_cr_r = V_cr_r(:,1);                     % select its eigenshape

%% Control reversal for the stabilizer
% clear fb, Lq, Lb of the rudder
stab = ttail;
stab.b(1).fb = zeros(size(stab.b(1).fb));
stab.b(1).Lq = zeros(size(stab.b(1).Lq));
stab.b(1).Lb = zeros(size(stab.b(1).Lb));

stab = m_compute_matrices(stab);

% Assebly stiffness the K matrices for cntr_rev
K_cr_s = [stab.K, zeros(size(stab.K,1),1); zeros(1,size(stab.K,2)),0];
Ka_cr_s = [stab.Ka, stab.fb; stab.Lq, stab.Lb];

% Find the control reversal (cr) dynamic pressure
[V_cr_s,D_cr_s]= eig(full(K_cr_s),full(Ka_cr_s));
q_cr_s = diag(D_cr_s);
q_cr_s = q_cr_s.*(abs(imag(q_cr_s))<10^-3);
[q_cr_s,I] = sort(real(q_cr_s));
V_cr_s = V_cr_s(:,I);                         % sort the eigenshapes
V_cr_s(:,q_cr_s<1)=[];                    % select the eigenshapes with positive eig
q_cr_s(q_cr_s<1)=[];
q_cr_s = q_cr_s(1);                   % select the minimum positive q_inf
V_cr_s = V_cr_s(:,1);                     % select its eigenshape


%% Plots
if 1
    % switch on the aero properties for the plot
    phi = rad2deg(0);
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    for g = 1:length(V_cr_s)-1
        if mod(g,6)==3 || mod(g,6)==5
            V_cr_s(g) = -V_cr_s(g);
        end
    end
    m_plot_eigenshape(ttail,options,((V_cr_s(1:end-1)))*2);
    
end

[T, a, P, rho] = atmosisa(0:100:11000);
v_cr = sqrt(q_cr_s*2./rho);
M_cr = v_cr./a;
if 1 % different figures
    figure(3)
    plot(0:100:11000,v_cr,'LineWidth',2)
    grid on
    xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
    ylabel('Velocity $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    set(gcf, 'Position',  [40, 40, 300, 300])
    saveas(figure(3),'divergence_tail_1','epsc')
    
    figure(4)
    plot(0:100:11000,M_cr,'LineWidth',2)
    grid on
    xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
    ylabel('Mach [-]','fontsize',14,'interpreter','latex')
    set(gcf, 'Position',  [40, 40, 300, 300])
    saveas(figure(4),'divergence_tail_2','epsc')
end

if 1
    % switch on the aero properties for the plot
    phi = rad2deg(0);
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    for g = 1:length(V_cr_r)-1
        if mod(g,6)==3 || mod(g,6)==5
            V_cr_r(g) = -V_cr_r(g);
        end
    end
    m_plot_eigenshape(ttail,options,((V_cr_r(1:end-1)))*3);
    
end

[T, a, P, rho] = atmosisa(0:100:11000);
v_cr = sqrt(q_cr_r*2./rho);
M_cr = v_cr./a;

if 1 % different figures
    figure(10)
    plot(0:100:11000,v_cr,'LineWidth',2)
    grid on
    xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
    ylabel('Velocity $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    set(gcf, 'Position',  [40, 40, 300, 300])
    saveas(figure(10),'cnt_rev_tail_rud_1','epsc')
    
    figure(11)
    plot(0:100:11000,M_cr,'LineWidth',2)
    grid on
    xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
    ylabel('Mach [-]','fontsize',14,'interpreter','latex')
    set(gcf, 'Position',  [40, 40, 300, 300])
    saveas(figure(11),'cnt_rev_tail_rud_2','epsc')
end

