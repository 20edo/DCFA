function b=b_init()
% This function initializes a beam
% (Geometry)
% b.cart                      % True if the section is defined as z=@(y)              [bool]
% b.Zmin                      % Coordinate of the 'lower' boundary as a function of x and y
% b.Zmax                      % Coordinate uf the 'upper boundary as a function of x and y
% b.Ymin                      % Coordinate of the 'lower' y boundary as a function of x
% b.Ymax                      % Coordinate of the 'upper' y boundary as a function of x
% b.pol                       % True if the section is defined in polar coordinates
% b.Rhomin                    % Coordinate of the 'lower' rho boundary as a function of x and th
% b.Rhomax                    % Coordinate of the 'upper' rho boundary as a function of x and th
% b.Thmin                     % Coordinate of the 'lower' th boundary as a function of x
% b.Thmax                     % Coordinate of the 'upper' th boundary as a function of x
% % All functions must be defined in the same reference of the geometry (pol
% % or cart)
% b.E                         % Young modulus of the point in the beam
% b.G                         % Shear modulus of the pooint in the beam
% b.rho                       % Density of the point of the beam
% % (Properties)
% b.L                         % Length of the beam     [m]
% b.in                        % Vector of the internal nodes of the beam
% b.nel                       % number of elements
% % (matrices)
% b.M                         % Mass matrix of the bar
% b.K                         % Stiffness matrix of the bar
% % Geometry of the beam in 3d space
% b.o                         % Coordinates of the origin of the beam              [3*1]
% b.v                         % Coordinates of the versor along the beam develops  [3*1]
% % Constraint are specified by a vector that contains true if the dof is
% % constrained and false if the dof is not.
% b.oc                        % Constraint at the origin                           [6*1]
% b.ec                        % Constaint at the end                               [6*1]
% b.on                        % Number associated to the node at the origin
% b.en                        % Number associated to the node at the end
% b.name                      % Optional name of the beam (string)
% functions defined on b:
% b_build_elemnts

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
%% 


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
b.v=nan(3,1);              % Coordinates of the versor along the beam develops  [3*1]
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
