clear all; clc; close all;
cd ..
cd ..
cd generate_model\
inside = 0;
pos_x = [-2:-0.5:-7]; 
pos_z = [-0.5:-0.5:-6]; 
pos = []; 
for i = 1:length(pos_x) 
    for j = 1:length(pos_z)
        pos = [pos; pos_x(i), pos_z(j)]; 
    end
end
pos = [pos; 1,1]; 
max = size(pos,1)-1;
save pos
cd ..
cd aero_analysis\Steady\
velocity = zeros(max+1,1);
velocity(end) = 1;
save velocity
%% Cycle in order to make:
% discharge_from_inside
% discharge_from_outside
if 0
    contator = 1;
    while contator<=length(velocity)-1
        flutter_study2
        flutter_velocity = v(i);
        load velocity
        contator = velocity(end)
        velocity(contator) = flutter_velocity;
        contator = contator + 1;
        velocity(end) = contator;
        clearvars -except velocity contator
        save velocity
    end
end
%% 
clear all 
load flutter_engines
pos_x = [-2:-0.5:-7]; 
pos_z = [-0.5:-0.5:-6]; 
velocity = reshape(velocity(1:end-1),length(pos_z),length(pos_x)); 
velocity = [pos_z',velocity]; 
velocity = [velocity; 0, pos_x]; 

