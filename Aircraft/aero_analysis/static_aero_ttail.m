% Static aero analysis of the T-Tail
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
clear all, close all, clc
cd ..
%% Generate the aircraft model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\


%% Build the wing model
ttail=m_init();
ttail.en=[en_ground(aircraft.en(5).x)];
ttail_list=[aircraft.b(13) aircraft.b(14) aircraft.b(15)];
for i=1:length(ttail_list)
    ttail=m_add_beam(ttail,ttail_list(i));
end

%% Add the aero loads
ttail = m_add_aero_loads(ttail,[1,0,0]');

%% Find the divergence dynamic pressure

[V_div,D_div]= eigs(ttail.K,ttail.Ka,30,'smallestabs');
q_div = diag(D_div);
[q_div,I] = sort(real(q_div));
V_div = V_div(:,I);                         % sort the eigenshapes
V_div(:,q_div<0)=[];                    % select the eigenshapes with positive eig
q_div(q_div<0)=[];
q_div = q_div(1);                   % select the minimum positive q_inf
V_div = V_div(:,1);                     % select its eigenshape

%% Calculations for the plotting VTAS and MACH when altitude changes
[T, a, P, rho] = atmosisa(0:100:11000);
v_div = sqrt(q_div*2./rho);
M_div = v_div./a;

%% Find the control reversal dynamic pressure 
K = [ttail.K, zeros(size(ttail.K,1),1); zeros(1,size(ttail.K,2)),0]; 
Ka = [ttail.Ka, ttail.fb; ttail.Lq, ttail.Lb]; 

%% Find the control reversal (cr) dynamic pressure
[V_cr,D_cr]= eig(full(K),full(Ka));
q_cr = diag(D_cr);
[q_cr,I] = sort(real(q_cr));
V_cr = V_cr(:,I);                         % sort the eigenshapes
V_cr(:,q_cr<0)=[];                    % select the eigenshapes with positive eig
q_cr(q_cr<0)=[];
q_cr = q_cr(1);                   % select the minimum positive q_inf
V_cr = V_cr(:,1);                     % select its eigenshape

%% Plote and save results
if 0
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;         % we have only one eig
    m_plot_eigenshape(ttail,options,1/4*V_div)
    
    figure(2)
    title('Static divergence T-Tail modeshape')
    fig = figure(2)
    for h = 1:4
        if h==1
            fname = ['Static_Divergence_T-Tail_view3D'];
            saveas(fig,fname,'svg')
        elseif h == 2
            view([1 0 0])
            fname = ['Static_Divergence_T-Tail_viewX'];
            saveas(fig,fname,'svg')
        elseif h==3
            view([0 1 0])
            fname = ['Static_Divergence_T-Tail_viewY'];
            saveas(fig,fname,'svg')
        elseif h==4
            view([0 0 1])
            fname = ['Static_Divergence_T-Tail_viewZ'];
            saveas(fig,fname,'svg')
        end
    end
    
    
    figure(3)
    set(gcf, 'Position',  [40, 40, 700, 500])
    subplot(1,2,1)
    plot(0:100:11000,v_div,'LineWidth',2)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence VTAS [m/s]')
    title('Divergence VTAS T-Tail')
    subplot(1,2,2)
    plot(0:100:11000,M_div,'LineWidth',2)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Mach [-]')
    title('Divergence Mach T-Tail')
    saveas(figure(3),'Static_divergence_T-Tail_graph','svg')
end