function sc=sc_constant_p_circle(R,E,G,rho)

% Returns a circular section of radius R and constant properties
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

sc=sc_init();
sc.cart=false;              % True if the section is defined as z=@(y)              [bool]
sc.pol=true;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) 0.*th;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) R+0.*th;        % Coordinate of the 'upper' rho boundary
sc.Thmin=0;               % Coordinate of the 'lower' th boundary
sc.Thmax=2*pi;               % Coordinate of the 'upper' th boundary


% All functions must be defined in the same reference of the geometry
sc.E=@(Rho,th) E+0.*Rho+0.*th;
sc.G=@(Rho,th) G+0.*Rho+0.*th;
sc.rho=@(Rho,th) rho+0.*Rho+0.*th;
sc=sc_compute_property(sc);

end
