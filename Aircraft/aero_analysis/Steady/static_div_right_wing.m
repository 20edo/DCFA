% Static aero analysis of the clamped wing
% - Divergence
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
clear all , close all, clc
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
% wing.en=[en_ground(aircraft.en(7).x)];
% wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9)];
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end
%% Generate the straight wing model
wing_straight = wing;
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');
wing = m_add_mass_forces(wing);
wing = m_compute_matrices(wing);
wing_straight = m_add_aero_loads_straight(wing_straight,[1,0,0]');

%% Find the divergence dynamic pressure - swept wing
% we take the smallest one positive
[V,D]= eig(full(wing.K),full(wing.Ka));
q_div = diag(D);
q_div = q_div.*(abs(imag(q_div))<10^-3);
[q_div,I] = sort(real(q_div));
V = V(:,I);                             % sort the eigenshapes
V(:,q_div<1)=[];                        % select the eigenshapes with positive eig
q_div(q_div<1)=[];
q_div = q_div(1);                       % select the minimum positive q_inf
V_div = V(:,1);                         % select its eigenshape

%% Find the divergence dynamic pressure - straight wing
[V_straight,D_straight]= eig(full(wing_straight.K),full(wing_straight.Ka));
q_div_straight = diag(D_straight);
[q_div_straight,I] = sort(real(q_div_straight));
V_straight = V_straight(:,I);                             % sort the eigenshapes
V_straight(:,q_div_straight<1)=[];                        % select the eigenshapes with positive eig
q_div_straight(q_div_straight<1)=[];
q_div_straight = q_div_straight(1);                       % select the minimum positive q_inf
V_div_straight = V_straight(:,1);

%% Assebly stiffness the K matrices for cntr_rev
K_cr = [wing_straight.K, zeros(size(wing_straight.K,1),1); zeros(1,size(wing_straight.K,2)),0];
Ka_cr = [wing_straight.Ka, wing_straight.fb; wing_straight.Lq, wing_straight.Lb];

%% Find the control reversal (cr) dynamic pressure
[V_cr,D_cr]= eig(full(K_cr),full(Ka_cr));
q_cr = diag(D_cr);
[q_cr,I] = sort(real(q_cr));
V_cr = V_cr(:,I);                         % sort the eigenshapes
V_cr(:,q_cr<1)=[];                    % select the eigenshapes with positive eig
q_cr(q_cr<1)=[];
q_cr = q_cr(1);                   % select the minimum positive q_inf
V_cr = V_cr(:,1);                     % select its eigenshape


%% Calculations for the plotting VTAS and MACH when altitude changes - swept wing
[T, a, P, rho] = atmosisa(0:100:11000);
v = sqrt(q_div*2./rho);
M = v./a;

%% Calculations for the plotting VTAS and MACH when altitude changes - straight wing
[T_straight, a_straight, P_straight, rho_straight] = atmosisa(0:100:11000);
v_straight = sqrt(q_div_straight*2./rho_straight);
M_straight = v_straight./a_straight;

%% Plot and save results
close all
if 1
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
        m_plot_eigenshape(wing,options,(-4*(V_div)));
    end
    saveas(figure(2),'divergence_wing_5','epsc')
    %     figure(2)
    %     title('Static divergence Wing modeshape')
    %     fig = figure(2);
    %     for h = 1:4
    %         if h==1
    %             fname = ['Static_Divergence_Wing_view3D'];
    %             saveas(fig,fname,'svg')
    %         elseif h == 2
    %             view([1 0 0])
    %             fname = ['Static_Divergence_Wing_viewX'];
    %             saveas(fig,fname,'svg')
    %         elseif h==3
    %             view([0 1 0])
    %             fname = ['Static_Divergence_Wing_viewY'];
    %             saveas(fig,fname,'svg')
    %         elseif h==4
    %             view([0 0 1])
    %             fname = ['Static_Divergence_Wing_viewZ'];
    %             saveas(fig,fname,'svg')
    %         end
    %     end
    if 0 % just one figure
        figure(3)
        set(gcf, 'Position',  [40, 40, 700, 500])
        subplot(2,2,1)
        plot(0:100:11000,v,'LineWidth',2)
        grid on
        xlabel('Altitude [m]')
        ylabel('Divergence Velocity [m/s]')
        title('Divergence Velocity Swept Wing')
        subplot(2,2,2)
        plot(0:100:11000,M,'LineWidth',2)
        grid on
        xlabel('Altitude [m]')
        ylabel('Divergence Mach [-]')
        title('Divergence Mach Swept Wing')
        subplot(2,2,3)
        plot(0:100:11000,v_straight,'LineWidth',2)
        grid on
        xlabel('Altitude [m]')
        ylabel('Divergence Velocity [m/s]')
        title('Divergence Velocity Straight Wing')
        subplot(2,2,4)
        plot(0:100:11000,M_straight,'LineWidth',2)
        grid on
        xlabel('Altitude [m]')
        ylabel('Divergence Mach [-]')
        title('Divergence Mach Straight Wing')
        saveas(figure(3),'Static_divergence_Wing_graph','svg')
    end
    if 1 % save in differen plots
        figure(4)
        plot(0:100:11000,v,'LineWidth',2)
        grid on
        xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
        ylabel('Velocity $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
        set(gcf, 'Position',  [40, 40, 300, 300])
        saveas(figure(4),'divergence_wing_1','epsc')
        
        figure(5)
        plot(0:100:11000,M,'LineWidth',2)
        grid on
        xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
        ylabel('Mach [-]','fontsize',14,'interpreter','latex')
        set(gcf, 'Position',  [40, 40, 300, 300])
        saveas(figure(5),'divergence_wing_2','epsc')
        
        figure(6)
        plot(0:100:11000,v_straight,'LineWidth',2)
        grid on
        xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
        ylabel('Velocity $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
        set(gcf, 'Position',  [40, 40, 300, 300])
        saveas(figure(6),'divergence_wing_3','epsc')
        
        figure(7)
        plot(0:100:11000,M_straight,'LineWidth',2)
        grid on
        xlabel('Altitude $\quad$ $[m]$','fontsize',14,'interpreter','latex')
        ylabel('Mach [-]','fontsize',14,'interpreter','latex')
        set(gcf, 'Position',  [40, 40, 300, 300])
        saveas(figure(7),'divergence_wing_4','epsc')
    end
end