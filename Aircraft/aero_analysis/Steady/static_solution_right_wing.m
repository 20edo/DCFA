% Static aero analysis of the clamped wing
% - Static solution
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

%% Direct problem
% study the deflection of the wing 
[T, a, P, rho] = atmosisa(10000);
M = 0.7; 
v = M*a; 
q = 1/2*rho*v.^2; 

A = wing.K - q*wing.Ka;
b = q*wing.fa + wing.f; 
q_stat = A\b; 

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
    options.plotAlpha              = 0.8;
    
m_plot_eigenshape(wing,options,q_stat);
end
% saveas(figure(2),'direct_solution_1','epsc')
%%
q_stat = [zeros(6,1); q_stat]; 
beam1 = wing.b(1).A'*q_stat;
beam2 = wing.b(2).A'*q_stat;
beam3 = wing.b(3).A'*q_stat;
flex = [beam1(3:6:end); beam2(9:6:end); beam3(9:6:end)]; 
tor = rad2deg([beam1(4:6:end); beam2(10:6:end); beam3(10:6:end)]);  
figure(5)
yyaxis left
plot(linspace(0,28,length(flex)),flex,'LineWidth',2);
ylabel('Displacements $\quad$ $[m]$','fontsize',14,'interpreter','latex')
yyaxis right 
plot(linspace(0,28,length(tor)),tor,'LineWidth',2);
ylabel('Rotations $\quad$ $[deg]$','fontsize',14,'interpreter','latex')
xlabel('Wingspan $\quad$ $[m]$','fontsize',14,'interpreter','latex')
legend('Vertical displacements','Torsional rotation','fontsize',10,'interpreter','latex','Location','northwest')
grid on 
saveas(figure(5),'direct_solution_2','epsc')
