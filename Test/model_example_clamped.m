% This is a commented example of structure built with this program

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%

% Add to path the necessary folders 

cd ..
init
cd Test\
close all, clear all, clc
%% Initialize the model
L_shaped_structure=m_init();

%% Set the external node and add it to the model
node_1=en_ground([0,0,0]);
L_shaped_structure.en=node_1;

%% Build and add the first beam to the model
beam_1=b_constant_p_square(1,0.1,70*1e9,27*1e9,2700,30);  % Build beam
beam_1.name='beam_1';       % Give a name to the beam (optional)
beam_1.o=[0,0,0]';          % Set the origin of the first beam
beam_1.vx=[1,0,0]';         % Set the versor of the first beam 
                            % (aligned with the global x axis)
beam_1.vy=[0,1,0]';         % Set the versor of the y local axis of the beam
                            
beam_1.oc=true(6,1);        % The first beam is clamped at the origin (wit\h node_1==ground)
beam_1.ec=true(6,1);        % The first beam is clamped at the end (with beam_2)
L_shaped_structure=m_add_beam(L_shaped_structure,beam_1);   % Add beam to the model
Straigth_structure = L_shaped_structure; 


%% Build and add the second beam to the model
beam_2=b_constant_p_square(2,0.1,70*1e9,27*1e9,2700,30);
beam_2.name='beam_2';
beam_2.o=beam_1.o+beam_1.L*beam_1.vx;   % Beam 2 origin is coincident with the end of beam 1
beam_2.vx=[0,1,0]';                     % Beam 2 is aligned with the global y axis
beam_2.vy=[0,0,1]';                     % Set the versor of the y local axis of the beam
beam_2.oc=true(6,1);                    % Beam 2 is clamped at the origin
beam_2.ec=false(6,1);                   % Beam 2 is free at the end
L_shaped_structure=m_add_beam(L_shaped_structure,beam_2);

%% Build a straight element 
beam_4=b_constant_p_square(2,1,70*1e9,27*1e9,2700,30);
beam_4.name='beam_2';
beam_4.o=beam_1.o+beam_1.L*beam_1.vx;   % Beam 2 origin is coincident with the end of beam 1
beam_4.vx=[1,0,0]';                     % Beam 2 is aligned with the global y axis
beam_4.vy=[0,1,0]';                     % Set the versor of the y local axis of the beam
beam_4.oc=true(6,1);                    % Beam 2 is clamped at the origin
beam_4.ec=false(6,1);                   % Beam 2 is free at the end
Straigth_structure=m_add_beam(Straigth_structure,beam_4);
%% Add loads

L_shaped_structure=m_add_mass_forces(L_shaped_structure);

%% Compute matrices

L_shaped_structure=m_compute_matrices(L_shaped_structure);

%% Solve static problem

L_shaped_structure=m_static_solution(L_shaped_structure);

%% Build a T structure with a lumped mass at the -y end

beam_3=b_constant_p_square(2,0.30,70*1e9,27*1e9,2700,30);
beam_3.name='beam_3';
beam_3.o=beam_1.o+beam_1.L*beam_1.vx;   % Beam 3 origin is coincident with the end of beam 1
beam_3.vx=[0,-1,0]';                    % Beam 3 is aligned with the global -y axis
beam_3.vy=[0,0,1]';                     % Set the versor of the y local axis of the beam
beam_3.oc=true(6,1);                    % Beam 3 is clamped at the origin
beam_3.ec=true(6,1);                    % Beam 3 is clamped (with a lumped mass) at the end


lumped_mass=en_mass(beam_3.o+beam_1.L*beam_3.vx,1e3);    % Create a lumped mass external node 1e3 Kg
T_shaped_structure=L_shaped_structure;
T_shaped_structure.en=[T_shaped_structure.en lumped_mass];

T_shaped_structure=m_add_beam(L_shaped_structure,beam_3);


%% Compute matrices

T_shaped_structure=m_compute_matrices(T_shaped_structure);

model = T_shaped_structure;

%% prova
% cd .. 
% cd Model
% 
% n = size(model.M,1); 
% y0 = zeros(2*n,1); 
% % f vorrei che avesse dimenzioni iniziali 
% deg = 6;
% f = @(t) [zeros(deg-1,1); 5*1e7*(t>0) ; zeros(size(model.M,1)-deg,1)]; 
% model = m_solution_dynamic_clamped(...
%     model,[-1,5000],y0,f,50);
% 
% 
% m_plot_easy(model)












% cd ..
% cd .. 
% cd Model
options.plot_original          = 1;
options.plot_deformed          = 1;
options.plotColor              = 'green';
options.saveSTL                = 0;
options.point_section          = 2;
[fig] = m_plot3d(Straigth_structure,options)
% hold on 
% quiver3(0,0,0,5,0,0)
% quiver3(0,0,0,0,5,0)
% quiver3(0,0,0,0,0,5)




