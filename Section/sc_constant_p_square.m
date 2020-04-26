function sc=sc_constant_p_square(l,E,G,rho)

% Returns a square section of side l and constant properties
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
sc.cart=true;
sc.Zmax=@(y) l/2+y.*0;
sc.Zmin=@(y) -l/2+y.*0;
sc.Ymax=l/2;
sc.Ymin=-l/2;
sc.pol=false;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) nan;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) nan;        % Coordinate of the 'upper' rho boundary
sc.Thmin=nan;               % Coordinate of the 'lower' th boundary
sc.Thmax=nan;               % Coordinate of the 'upper' th boundary
% All functions must be defined in the same reference of the geometry
sc.E=@(y,z) E+0.*y+0.*z;
sc.G=@(y,z) G+0.*y+0.*z;
sc.rho=@(y,z) rho+0.*y+0.*z;
sc=sc_inertia(sc);


end
