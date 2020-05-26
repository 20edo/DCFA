% This script builds the model of the aircraft
% The origin is set at the nose of the aircraft, the x axis is aligned with
% the fuselage and points to the front

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

% payload

wings

% fuel

Ttail

engines

%% Compute properties

aircraft=m_compute_matrices(aircraft);

%% Plot

if 0
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    [fig] = m_plot3d(aircraft,options);
end


