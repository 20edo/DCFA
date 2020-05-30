function el=el_build_constant_p_square(b,insx,indx)
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
sc=sc_init();
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

%% Define constant section variables
l=sc.Ymax-sc.Ymin;
E=sc.E(0,0);
G=sc.G(0,0);
rho=sc.rho(0,0);
%% Compute section properties
sc=sc_constant_p_square(l,E,G,rho);
%% Set element properties
el = el_init();
el.sc=sc;
el.L=abs(indx.x-insx.x);
% el.sc.Iy = 0; 
% el.sc.Iz = 0; 
el=el_mass_assembly(el);
el=el_stiff_assembly(el);
el.C=zeros(size(el.M));
end
