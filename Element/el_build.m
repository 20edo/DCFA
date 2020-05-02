function el=el_build(b,insx,indx)
% Builds the element between two nodes, given the beam property.
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
%% Find the coordinate in which the section must be calculated
x=(insx.x+indx.x)/2;


%% Initialise section and set geometric and material properties

sc.cart=b.cart;              % True if the section is defined as z=@(y)              [bool]
sc.Zmin=@(y) b.Zmin(x,y);           % Coordinate of the 'lower' boundary
sc.Zmax=@(y) b.Zmax(x,y);           % Coordinate uf the 'upper boundary
sc.Ymin= b.Ymin(x);               % Coordinate of the 'lower' x boundary
sc.Ymax= b.Ymax(x);               % Coordinate of the 'upper' x boundary
sc.pol=b.pol;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) b.Rhomin(x,th);        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) b.Rhomax(x,th);        % Coordinate of the 'upper' rho boundary
sc.Thmin=b.Thmin(x);               % Coordinate of the 'lower' th boundary
sc.Thmax=b.Thmax(x);               % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
sc.E=@(y,z) b.E(x,y,z);            % Young modulus of the point in the section
sc.G=@(y,z) b.G(x,y,z);            % Shear modulus of the pooint in the section
sc.rho=@(y,z) b.rho(x,y,z);          % Density of the point of the section
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
sc.Jp=nan;                  % Polar moment of inertia [kg*mm^2]


%% Compute section properties
sc=sc_compute_property(sc);
sc=sc_torsion(sc);

%% Set element properties
el.sc=sc;
el.L=abs(indx.x-insx.x);
el.M=el_mass_assembly(el);
el.K=el_stiff_assembly(el);
end


