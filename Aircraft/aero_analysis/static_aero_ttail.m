% Static aero analysis of the T-Tail
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

ttail=m_init();
ttail.en=[en_ground(aircraft.en(5).x)];
ttail_list=[aircraft.b(13) aircraft.b(14) aircraft.b(15)];
  
for i=1:length(ttail_list)
    ttail=m_add_beam(ttail,ttail_list(i));
end

% Add the aero loads
ttail = m_add_aero_loads(ttail,[1,0,0]');

% Find the divergence dynamic pressure
q_div = eigs(ttail.K,ttail.Ka,30,'smallestabs');
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