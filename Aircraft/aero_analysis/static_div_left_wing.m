% Static aero analysis of the clamped wing
% - Divergence
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
% Move to the analysis folder
cd aero_analysis\

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
%% Generate the straight wing model
wing_straight = wing;
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');
wing_straight = m_add_aero_loads_straight(wing_straight,[1,0,0]');

%% Find the divergence dynamic pressure - swept wing
[V,D]= eig(full(wing.K),full(wing.Ka));
q_div = diag(D);
[q_div,I] = sort(real(q_div));
V = V(:,I);                             % sort the eigenshapes
V(:,q_div<0)=[];                        % select the eigenshapes with positive eig
q_div(q_div<0)=[];
q_div = q_div(1);                       % select the minimum positive q_inf
V_div = V(:,1);                         % select its eigenshape

%% Find the divergence dynamic pressure - straight wing
q_div_straight = eigs(wing_straight.K,wing_straight.Ka,30,'smallestabs');
q_div_straight = sort(real(q_div_straight));
q_div_straight(q_div_straight<0)=[];
q_div_straight = q_div_straight(1);


%% Calculations for the plotting VTAS and MACH when altitude changes - swept wing
[T, a, P, rho] = atmosisa(0:100:11000);
v = sqrt(q_div*2./rho);
M = v./a;

%% Calculations for the plotting VTAS and MACH when altitude changes - straight wing
[T_straight, a_straight, P_straight, rho_straight] = atmosisa(0:100:11000);
v_straight = sqrt(q_div_straight*2./rho_straight);
M_straight = v_straight./a_straight;

%% Plot and save results
if 1
    
    % switch on the aero properties for the plot 
    for i=4:5
        wing.b(i).ssh = true; 
    end
    
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    m_plot_eigenshape(wing,options,V_div/20)
    
    
    figure(2)
    title('Static divergence Wing modeshape')
    fig = figure(2)
    for h = 1:4
        if h==1
            fname = ['Static_Divergence_Wing_view3D'];
            saveas(fig,fname,'svg')
        elseif h == 2
            view([1 0 0])
            fname = ['Static_Divergence_Wing_viewX'];
            saveas(fig,fname,'svg')
        elseif h==3
            view([0 1 0])
            fname = ['Static_Divergence_Wing_viewY'];
            saveas(fig,fname,'svg')
        elseif h==4
            view([0 0 1])
            fname = ['Static_Divergence_Wing_viewZ'];
            saveas(fig,fname,'svg')
        end
    end
    
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