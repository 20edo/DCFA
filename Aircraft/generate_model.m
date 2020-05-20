% This script builds the model of the aircraft
% The origin is set at the nose of the aircraft, the x axis is align with
% the fuselage ond points to the front

%%
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


close all, clear all, clc

%% Add to path what needed

cd ..

addpath('Section')
addpath('Element')
addpath('Beam')
addpath('internal_nodes')
addpath('External_nodes')
addpath('Model')

cd Aircraft\

%% Define the nodes of the beam
Node1=en_free([0,0,0]);     % Where the node begins
Node2=en_free([-15 0 0]);    % Where the section from conical becomes cylindrical
Node3=en_free([-50 0 0]);    % Where the wings are clamped
Node4=en_free([-100 0 0]);   % Where the conical section of the tail begins
Node5=en_free([-120 0 3]); % Where the rudder is clamped
Node6=en_free([-140 0 6.5]);   % Where the aircraft ends

%% Initialize model

aircraft=m_init();
aircraft.en=[Node1, Node2, Node3, Node4, Node5, Node6];

% clear unuseful variables
clear Node1
clear Node2
clear Node3
clear Node4
clear Node5
clear Node6
clear i

%% Build the beams
fuselage

for i=1:length(fuselage_beams)
    aircraft=m_add_beam(aircraft,fuselage_beams(i));
end

%% Compute properties

aircraft=m_compute_matrices(aircraft);

options.plot_original          = 1;
options.plot_deformed          = 1;
options.plotColor              = 'green';
options.saveSTL                = 0;
options.point_section          = 4;
[fig] = m_plot3d(aircraft,options)




