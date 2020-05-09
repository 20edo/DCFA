function sc=sc_constant_p_square(l,E,G,rho)

% Returns a square section of side l and constant properties.
% Torsion properties are included.
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
sc.cart=true;
sc.Zmax=@(y) l/2+y.*0;
sc.Zmin=@(y) -l/2+y.*0;
sc.Ymax=l/2;
sc.Ymin=-l/2;
sc.pol=false;               % True if the section is defined in polar coordinates
% All functions must be defined in the same reference of the geometry
sc.E=@(y,z) E+0.*y+0.*z;
sc.G=@(y,z) G+0.*y+0.*z;
sc.rho=@(y,z) rho+0.*y+0.*z;
%% Inertia
sc.m=rho*l.^2;
sc.ycg=0;
sc.zcg=0;
sc.Iy=l^4/12;
sc.Iz=l^4/12;
sc.Iyz=0;
sc.Jp=sc.m*l^4/6;
%% Elastic
sc.GA=G.*l^2;
sc.EA=E.*l^2;
sc.ya=0;
sc.za=0;
sc.EJy=E*sc.Iy;
sc.EJz=E*sc.Iz;
sc.EJyz=0;
sc.GJ=2.25*(l/2)^4*G;
sc.yct=0;
sc.zct=0;
end
