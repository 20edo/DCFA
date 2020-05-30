% Static aero analysis of the clamped aircraft
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

% Add the aero loads
aircraft = m_add_aero_loads(aircraft,[1,0,0]');

% Find the divergence dynamic pressure
q_div = eigs(aircraft.K(7:end,7:end),aircraft.Ka(7:end,7:end),30,'smallestabs');
q_div = sort(real(q_div));
q_div(q_div<0)=[]; 
q_div = q_div(1)