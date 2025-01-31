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



%% Compute section properties
% First of all, properties wrt (0,0) are calculated to obtain the correct
% yct and zct, then the section is moved and properties are calculated wrt
% ct, finally yct and zct voordinates are stored in the appropriate field
sc=sc_compute_property(sc); % Calculate inertia and elastic properties
sc=sc_torsion(sc);          % Calculate yct and zct
sc_ct=sc_move_geo_reference(sc); 
sc_ct=sc_compute_property(sc_ct);   % Calculate inertia and elastic properties
sc_ct=sc_torsion(sc_ct);    % Calculate GJ wrt ct
sc_ct.yct=sc.yct;
sc_ct.zct=sc.zct;



%% Set element properties
el.sc=sc_ct;
el.L=abs(indx.x-insx.x);
el=el_mass_assembly(el);
el=el_stiff_assembly(el);
el.C=zeros(size(el.M));
end


