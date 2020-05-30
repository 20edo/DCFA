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

% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');

% Find the divergence dynamic pressure
q_div = eigs(wing.K,wing.Ka,1,'smallestabs');