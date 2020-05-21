% This script builds the model of the aircraft
% The origin is set at the nose of the aircraft, the x axis is align with
% the fuselage ond points to the front

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


%% Initialize model

aircraft=m_init();

%% Build the beams

fuselage 

wings

Ttail

%% Compute properties

aircraft=m_compute_matrices(aircraft);

%% Plot


options.plot_original          = 1;
options.plot_deformed          = 1;
options.plotColor              = 'green';
options.saveSTL                = 0;
options.point_section          = 8;
[fig] = m_plot3d(aircraft,options)




