% This is a commented example of structure built with this program

% DCFA swept wing assignement
%
% Teamwork
% Team members: Pasturenzi Lorenzo    944610
%               Tacchi Alberto        944579
%               Venti Edoardo         944421
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
beam_1=b_constant_p_square(1e3,100,70*1e6,27*1e6,2700,30);  % Build beam
beam_1.name='beam_1';       % Give a name to the beam (optional)
beam_1.o=[0,0,0]';          % Set the origin of the first beam
beam_1.vx=[1,0,0]';         % Set the versor of the first beam 
                            % (aligned with the global x axis)
beam_1.vy=[0,0,1]';         % Set the versor of the y local axis of the beam
                            
beam_1.oc=true(6,1);        % The first beam is clamped at the origin (with node_1==ground)
beam_1.ec=true(6,1);        % The first beam is clamped at the end (with beam_2)
L_shaped_structure=m_add_beam(L_shaped_structure,beam_1);   % Add beam to the model

%% Build and add the second beam to the model
beam_2=b_constant_p_square(5e3,30,70*1e6,27*1e6,2700,30);
beam_2.name='beam_2';
beam_2.o=beam_1.o+beam_1.L*beam_1.v;    % Beam 2 origin is coincident with the end of beam 1
beam_2.vx=[0,1,0]';                     % Beam 2 is aligned with the global y axis
beam_2.vy=[0,0,1]';                     % Set the versor of the y local axis of the beam
beam_2.oc=true(6,1);                    % Beam 2 is clamped at the origin
beam_2.ec=false(6,1);                   % Beam 2 is free at the end
L_shaped_structure=m_add_beam(L_shaped_structure,beam_2);

%% Compute matrices

% L_shaped_structure=m_compute_matrices(L_shaped_structure);

%% Build a T structure with a lumped mass at the -y end

beam_3=b_constant_p_square(5e3,30,70*1e6,27*1e6,2700,30);
beam_3.name='beam_3';
beam_3.o=beam_1.o+beam_1.L*beam_1.v;    % Beam 3 origin is coincident with the end of beam 1
beam_3.vx=[0,-1,0]';                    % Beam 3 is aligned with the global -y axis
beam_3.vy=[0,0,1]';                     % Set the versor of the y local axis of the beam
beam_3.oc=true(6,1);                    % Beam 3 is clamped at the origin
beam_3.ec=true(6,1);                    % Beam 3 is clamped (with a lumped mass) at the end


lumped_mass=en_mass(beam_3.o+beam_1.L*beam_3.v,1e3);    % Create a lumped mass external node 1e3 Kg
T_shaped_structure=L_shaped_structure;
T_shaped_structure.en=[T_shaped_structure.en lumped_mass];

T_shaped_structure=m_add_beam(L_shaped_structure,beam_3);

%% Compute matrices

% T_shaped_structure=m_compute_matrices(L_shaped_structure);





