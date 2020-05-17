function b=b_constant_p_tube(L,R,t,E,G,rho,nel)
% Creates a beam with constant properties and tube section of side external radius R(x) and thickness t(x). 
% The length of the beam is L. nin is the number of internal nodes.
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
b=b_init();
% (Geometry)
b.cart=false;              % True if the section is defined as z=@(y)              [bool]
b.Zmin=@(x,y) nan;         % Coordinate of the 'lower' boundary
b.Zmax=@(x,y) nan;       % Coordinate uf the 'upper boundary
b.Ymin=@(x) nan;    % Coordinate of the 'lower' x boundary
b.Ymax=@(x) nan;   % Coordinate of the 'upper' x boundary
b.pol=true;               % True if the section is defined in polar coordinates
b.Rhomin=@(x,th) R(x)-t(x)+th.*0;      % Coordinate of the 'lower' rho boundary
b.Rhomax=@(x,th) R(x)+th.*0;      % Coordinate of the 'upper' rho boundary
b.Thmin=@(x) 0;          % Coordinate of the 'lower' th boundary
b.Thmax=@(x) 2*pi;          % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry (pol
% or cart)
b.E=@(x,y,z) E+0.*x+0.*y+0.*z;        % Young modulus of the point in the beam
b.G=@(x,y,z) G+0.*x+0.*y+0.*z;          % Shear modulus of the pooint in the beam
b.rho=@(x,y,z) rho+0.*x+0.*y+0.*z;       % Density of the point of the beam
% (Properties)
b.L=L;                   % Length of the beam     [m]
b.nel=nel;               % Number of elements

%% Speed up since the section is constant
b=b_build_elements_constant_p_tube(b);
b=b_build_matrices(b);
end
