function sc=sc_constant_p_tube(R,t,E,G,rho)

% Returns a circular tube of radius R and thickness t and constant properties
% Improve this function: write the exact values instad of integrating
% numerically
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
%% Data structures definition:

if t>=R
    error('The tube does not exist')
end
% section   
% Upper-case letter are referred to the geometry
% (Check ym in DCFA 4.26 is y cg)
sc.cart=false;              % True if the section is defined as z=@(y)              [bool]
sc.Zmin=nan;           % Coordinate of the 'lower' boundary
sc.Zmax=nan;           % Coordinate uf the 'upper boundary
sc.Ymin= nan;               % Coordinate of the 'lower' x boundary
sc.Ymax= nan;               % Coordinate of the 'upper' x boundary
sc.pol=true;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) (R-t)+0.*th;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) R+0.*th;        % Coordinate of the 'upper' rho boundary
sc.Thmin=0;               % Coordinate of the 'lower' th boundary
sc.Thmax=2*pi;               % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
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
sc.Jp=nan;                  % Polar moment of inertia
sc.cart=false;
sc.Zmax=nan;
sc.Zmin=nan;
sc.Ymax=nan;
sc.Ymin=nan;

% All functions must be defined in the same reference of the geometry
sc.E=@(Rho,th) E+0.*Rho+0.*th;
sc.G=@(Rho,th) G+0.*Rho+0.*th;
sc.rho=@(Rho,th) rho+0.*Rho+0.*th;
sc=sc_compute_property(sc);

end
