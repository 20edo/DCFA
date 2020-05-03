function b=b_constant_p_square_test(L,l,E,G,rho,nel)
% Creates a beam with constant properties and square section of side l. The
% length of the beam is L. nin is the number of internal nodes. It
% calculates a the laplacian for all section even though they are all
% equal. Use this function to test the improvements made in solving
% sections.
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
% Beam
% (Geometry)
b.cart=true;              % True if the section is defined as z=@(y)              [bool]
b.Zmin=@(x,y) -l/2+0.*x+0.*y;         % Coordinate of the 'lower' boundary
b.Zmax=@(x,y) l/2+0.*x+0.*y;         % Coordinate uf the 'upper boundary
b.Ymin= @(x) -l/2+0.*x;    % Coordinate of the 'lower' x boundary
b.Ymax= @(x) l/2+0.*x;    % Coordinate of the 'upper' x boundary
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
b.in=nan;                  % Vector of the internal nodes of the beam
b.nel=nel;               % Number of elements


b=b_build_elements(b);
b=b_build_matrices(b);
end
