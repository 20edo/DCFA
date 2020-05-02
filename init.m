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
el.mass=nan;    % Mass per unit length          [Kg/mm]
el.w=nan;       % Weigth per unit length        [N/mm]
el.EJy=nan;     % Elastic stiffness wrt y axis  [N*mm^2]
el.EJz=nan;     % Elastic stiffness wrt z axis  [N*mm^2]
el.EJyz=nan;    % Coupling elastic stiffness    [N*mm^2]
el.EA=nan;      % Elastic stiffness (traction)  [N]
el.GA=nan;      % Elastic stiffness (shear)     [N]
el.za=nan;      % z coordinate of elastic axis  [mm]
el.ya=nan;      % y coordinate of elastic axis  [mm]
el.zm=nan;      % z coordinate of inertia axis  [mm]
el.ym=nan;      % y coordinate of inertia axis  [mm]
el.Iy=nan;      % Inertia property wrt y axis   [Kg/mm]
el.Iz=nan;      % Inertia property wrt z axis   [Kg/mm]
el.Iyz=nan;     % Coupling inertia property     [Kg/mm]
el.GJ=nan;      % Torsional stiffness           [N*mm/rad]
el.Jp=nan;      % Polar inertia moment          [Kg*mm^4] (to check)
el.M=NaN(6);    % Mass matrix                   6*6 martix
el.K=NaN(6);    % Stiffness matrix              6*6 matrix
el.C=NaN(6);    % Dissipation matrix            6*6 matrix
% functions defined on el: 
%               

% Beam (to do)


% node
nd.init_pos=NaN(6,1);       % Initial position          [x,y,z,thx,thy,thz]
nd.displacement=NaN(6,1);   % Displacement vector       [sx,sy,sz,dthx,dthy,dthz]
nd.M=NaN(6);                % Lumped mass matrix        6*6 matrix
nd.K=NaN(6);                % Lumped stiffness matrix   6*6 matrix
nd.C=NaN(6);                % Lumped dissipation matrix 6*6 matrix
% functions defined on nd:



              