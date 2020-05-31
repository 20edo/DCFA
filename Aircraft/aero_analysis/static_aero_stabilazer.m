% Static aero analysis of the right stabilazer
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
% generate the wing model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\


% Build the wing model

stab=m_init();
stab.en=[en_ground(aircraft.en(14).x)];
stab_list=[aircraft.b(14)];
  
for i=1:length(stab_list)
    stab=m_add_beam(stab,stab_list(i));
end

% Add the aero loads
stab = m_add_aero_loads(stab,[1,0,0]');

% Find the divergence dynamic pressure
q_div = eigs(stab.K,stab.Ka,30,'smallestabs');
q_div = sort(real(q_div));
q_div(q_div<0)=[]; 
q_div = q_div(1)

[T, a, P, rho] = atmosisa(0:100:11000); 
v = sqrt(q_div*2./rho); 
M = v./a; 

subplot(1,2,1)
    plot(0:100:11000,v)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Velocity [m/s]')
    title('Divergence Velocity - Altitude')
subplot(1,2,2)
    plot(0:100:11000,M)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Mach [-]')
    title('Divergence Mach - Altitude')