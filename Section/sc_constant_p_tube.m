function sc=sc_constant_p_tube(R,t,E,G,rho)

% Returns a circular tube of radius R and thickness t and constant properties

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
%% Data structures definition:

if t>=R
    error('The tube does not exist')
end
% section 
sc=sc_init();
sc.pol=true;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) (R-t)+0.*th;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) R+0.*th;        % Coordinate of the 'upper' rho boundary
sc.Thmin=0;               % Coordinate of the 'lower' th boundary
sc.Thmax=2*pi;               % Coordinate of the 'upper' th boundary
sc.cart=false;
sc.E=@(Rho,th) E+0.*Rho+0.*th;
sc.G=@(Rho,th) G+0.*Rho+0.*th;
sc.rho=@(Rho,th) rho+0.*Rho+0.*th;
% sc=sc_compute_property(sc);

% Inertia properties
sc.ycg=0;
sc.zcg=0;
sc.Iy=pi/4*(R^4-(R-t)^4)/4;
sc.Iz=pi/4*(R^4-(R-t)^4)/4;
sc.Iyz=0;
sc.Jp=pi/2*(R^4-(R-t)^4)/4;
sc.m=(R^2-(R-t)^2)*pi*rho;

% Elastic properties
sc.GA=sc.G(0,0)*pi*(R^2-(R-t)^2);
sc.EA=sc.E(0,0)*pi*(R^2-(R-t)^2);
sc.ya=0;
sc.za=0;
sc.EJy=sc.E(0,0)*sc.Iy;
sc.EJz=sc.EJy;
sc.EJyz=0;

% Torsion properties
sc.yct=0;
sc.zct=0;
sc.GJ=sc.Jp*sc.G(0,0);

end
