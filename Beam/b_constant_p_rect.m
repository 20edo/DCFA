function b=b_constant_p_rect(L,l1,l2,E,G,rho,nel)
% Creates a beam with constant properties and rectangular section of side l1 along y and l2 along z.
% The length of the beam is L. nel is the number of elements.
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
% Beam
b=b_init();
% (Geometry)
b.cart=true;              % True if the section is defined as z=@(y)              [bool]
b.Zmin=@(x,y) -l2/2+0.*x+0.*y;         % Coordinate of the 'lower' boundary
b.Zmax=@(x,y) l2/2+0.*x+0.*y;         % Coordinate uf the 'upper boundary
b.Ymin= @(x) -l1/2+0.*x;    % Coordinate of the 'lower' x boundary
b.Ymax= @(x) l1/2+0.*x;    % Coordinate of the 'upper' x boundary
b.pol=false;               % True if the section is defined in polar coordinates
b.Rhomin=@(x,th) nan;      % Coordinate of the 'lower' rho boundary
b.Rhomax=@(x,th) nan;      % Coordinate of the 'upper' rho boundary
b.Thmin=@(x) nan;          % Coordinate of the 'lower' th boundary
b.Thmax=@(x) nan;          % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
b.E=@(x,y,z) E+0.*x+0.*y+0.*z;        % Young modulus of the point in the beam
b.G=@(x,y,z) G+0.*x+0.*y+0.*z;          % Shear modulus of the pooint in the beam
b.rho=@(x,y,z) rho+0.*x+0.*y+0.*z;        % Density of the point of the beam
% (Properties)
b.L=L;                   % Length of the beam     [m]
b.nel=nel;               % Number of elements

%% Speed up since the section is constant
b=b_build_elements_constant_p_square(b);
b=b_build_matrices(b);
end
