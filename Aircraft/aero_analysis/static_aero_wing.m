% Static aero analysis of the clamped wing
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

wing=m_init();
wing.en=[en_ground(aircraft.en(7).x) ...
    aircraft.en(17) aircraft.en(18)];
wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9) aircraft.b(16) aircraft.b(17)];
  
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end
wing_straight = wing; 
% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');
wing_straight = m_add_aero_loads_straight(wing_straight,[1,0,0]'); 

% Find the divergence dynamic pressure
q_div = eigs(wing.K,wing.Ka,30,'smallestabs');
q_div = sort(real(q_div));
q_div(q_div<0)=[]; 
q_div = q_div(1);

q_div_straight = eigs(wing_straight.K,wing_straight.Ka,30,'smallestabs');
q_div_straight = sort(real(q_div_straight));
q_div_straight(q_div_straight<0)=[]; 
q_div_straight = q_div_straight(1);



[T, a, P, rho] = atmosisa(0:100:11000); 
v = sqrt(q_div*2./rho); 
M = v./a; 

[T_straight, a_straight, P_straight, rho_straight] = atmosisa(0:100:11000); 
v_straight = sqrt(q_div_straight*2./rho_straight); 
M_straight = v_straight./a_straight; 

subplot(2,2,1)
    plot(0:100:11000,v)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Velocity [m/s]')
    title('Divergence Velocity Swept Wing')
subplot(2,2,2)
    plot(0:100:11000,M)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Mach [-]')
    title('Divergence Mach Swept Wing')
subplot(2,2,3)
    plot(0:100:11000,v_straight)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Velocity [m/s]')
    title('Divergence Velocity Straight Wing')
subplot(2,2,4)
    plot(0:100:11000,M_straight)
    grid on
    xlabel('Altitude [m]')
    ylabel('Divergence Mach [-]')
    title('Divergence Mach Straight Wing')




