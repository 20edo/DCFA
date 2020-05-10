% In this script constant and data structures are definied
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
%% Adding required repositories
addpath('Section')
addpath('Element')
addpath('Beam')
addpath('internal_nodes')
addpath('External_nodes')
addpath('Model')
%% Engineering constant definition
global g
g=9.81;        % Gravity acceleration          [m/s^2]

%% Data structures definition:

% section   
% Upper-case letter are referred to the geometry
% (Check ym in DCFA 4.26 is y cg)
sc.cart=false;              % True if the section is defined as z=@(y)              [bool]
sc.Zmin=@(y) nan;           % Coordinate of the 'lower' boundary
sc.Zmax=@(y) nan;           % Coordinate uf the 'upper boundary
sc.Ymin= nan;               % Coordinate of the 'lower' x boundary
sc.Ymax= nan;               % Coordinate of the 'upper' x boundary
sc.pol=false;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) nan;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) nan;        % Coordinate of the 'upper' rho boundary
sc.Thmin=nan;               % Coordinate of the 'lower' th boundary
sc.Thmax=nan;               % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
sc.E=@(y,z) nan;            % Young modulus of the point in the section
sc.G=@(y,z) nan;            % Shear modulus of the pooint in the section
sc.rho=@(y,z) nan;          % Density of the point of the section
sc.m= nan;                  % Mass of the section
sc.ycg=nan;                 % ycg of the section
sc.zcg=nan;                 % zcg of the section
sc.Iy=nan;                  % Inertia wrt y axis
sc.Iz=nan;                  % Inertia wrt z axis
sc.Iyz=nan;                 % Inertia (coupling term)
sc.EJy=nan;                 % Elastic stiffness wrt y axis  [N*mm^2]
sc.EJz=nan;                 % Elastic stiffness wrt z axis  [N*mm^2]
sc.EJyz=nan;                % Coupling elastic stiffness    [N*mm^2]
sc.EA=nan;                  % Elastic stiffness (traction)  [N]
sc.GA=nan;                  % Elastic stiffness (shear)     [N]
sc.za=nan;                  % z coordinate of elastic axis  [mm]
sc.ya=nan;                  % y coordinate of elastic axis  [mm]
sc.Iy=nan;                  % Inertia property wrt y axis   [Kg/mm]
sc.Iz=nan;                  % Inertia property wrt z axis   [Kg/mm]
sc.Iyz=nan;                 % Coupling inertia property     [Kg/mm]
sc.GJ=nan;                  % Torsional stiffness
sc.yct=nan;                 % Coordinate of the shear center[mm]
sc.zct=nan;                 % Coordinate of the shear center[mm] 
sc.Jp=nan;                  % Polar moment of inertia       [kg*mm^2]
% functions defined on sc:
% sc_inertia
% sc_elastic
% sc_compute_property
% sc_constant_p_square
% sc_torsion


% element

el.sc=sc;       % Section of the element
el.L=nan;       % Length                        [mm]
el.M=NaN(6);    % Mass matrix                   6*6 martix
el.K=NaN(6);    % Stiffness matrix              6*6 matrix
el.C=NaN(6);    % Dissipation matrix            6*6 matrix
% functions defined on el: 
% 

% Beam
% (Geometry)
b.cart=false;              % True if the section is defined as z=@(y)              [bool]
b.Zmin=@(x,y) nan;         % Coordinate of the 'lower' boundary
b.Zmax=@(x,y) nan;         % Coordinate uf the 'upper boundary
b.Ymin= @(x) nan;          % Coordinate of the 'lower' x boundary
b.Ymax= @(x) nan;          % Coordinate of the 'upper' x boundary
b.pol=false;               % True if the section is defined in polar coordinates
b.Rhomin=@(x,th) nan;      % Coordinate of the 'lower' rho boundary
b.Rhomax=@(x,th) nan;      % Coordinate of the 'upper' rho boundary
b.Thmin=@(x) nan;          % Coordinate of the 'lower' th boundary
b.Thmax=@(x) nan;          % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
b.E=@(x,y,z) nan;          % Young modulus of the point in the beam
b.G=@(x,y,z) nan;          % Shear modulus of the pooint in the beam
b.rho=@(x,y,z) nan;        % Density of the point of the beam
% (Properties)
b.L=nan;                   % Length of the beam     [m]
b.in=nan;                  % Vector of the internal nodes of the beam
b.nel=nan;                     % number of elements
% (matrices)
b.M=nan;                   % Mass matrix of the bar
b.K=nan;                   % Stiffness matrix of the bar
% Geometry of the beam in 3d space
b.o=nan(3,1);              % Coordinates of the origin of the beam              [3*1]
b.vx=nan(3,1);             % Coordinates of the versor along the beam develops  [3*1]
b.vy=nan(3,1);             % Coordinates of the y versor of the section         [3*1] 
% The following fields exist if and only if the beam is inside a model
% Constraint are specified by a vector that contains true if the dof is
% constrained and false if the dof is not.
b.oc=nan(6,1);             % Constraint at the origin                           [6*1]
b.ec=nan(6,1);             % Constaint at the end                               [6*1]
b.on=nan;                  % Number associated to the node at the origin
b.en=nan;                  % Number associated to the node at the end
b.name=[];                 % Optional name of the beam (string)
% functions defined on b:
% b_build_elemnts


% Internal node
in.x=nan;               % Position along the beam axis
in.d=zeros(6,1);        % Displacement vector [6*1] vector  [mm]

% External node
en.x=nan(3,1);              % Position in the 3d space [3*1]    [mm]
en.d=zeros(6,1);            % Displacement vector [6*1] vector  [mm]
en.c=nan(6,1);              % Constraint 
en.M=nan(6);                % Mass matrix
en.K=nan(6);                % Stiffness matrix
en.C=nan(6);                % Dissipation matrix


% model
m.en=[];                    % Extrenal nodes (where more than one beam join or
                            % masses, forces, constraints,... are applied)
m.b=[];
m.M=nan;                    % Mass matrix of the model
m.K=nan;                    % Stiffness matrix of the model


% % % % % node
% % % % nd.init_pos=NaN(6,1);       % Initial position          [x,y,z,thx,thy,thz]
% % % % nd.displacement=NaN(6,1);   % Displacement vector       [sx,sy,sz,dthx,dthy,dthz]
% % % % nd.M=NaN(6);                % Lumped mass matrix        6*6 matrix
% % % % nd.K=NaN(6);                % Lumped stiffness matrix   6*6 matrix
% % % % nd.C=NaN(6);                % Lumped dissipation matrix 6*6 matrix
% % % % % functions defined on nd:
% % % % 


              