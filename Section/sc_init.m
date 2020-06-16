function sc=sc_init()
% Initialize a section
% Functions defined on section object:
% sc_inertia
% sc_elastic
% sc_compute_property
% sc_constant_p_square
% sc_torsion



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


sc.navier=zeros(6);         % Has no sense for continuous sections, but need it to assemble the model
