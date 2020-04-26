function [sc]=sc_elastic(sc)

% Calculate the elastic propreties of a section
% The function handles must be defined to accept arrays and return arrays
% with the same dimension (trick: add y.*0+z.*0)
% Missing torsion
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

%% Cartesian coordinates
if sc.cart==true
    sc.GA=integral2(sc.G,sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.EA=integral2(sc.E,sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.ya=1/sc.EA*integral2(@(y,z) y.*sc.E(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.za=1/sc.EA*integral2(@(y,z) z.*sc.E(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.EJy=integral2(@(y,z) z.^2.*sc.E(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.EJz=integral2(@(y,z) y.^2.*sc.E(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.EJyz=integral2(@(y,z) y.*z.*sc.E(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    
%% Polar coordinates        y=r.*cos(t)/z=r.sin(t)
elseif sc.polar==true
    sc.GA=intergal2(@(r,t) sc.G(r,t).*r, sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.EA=intergal2(@(r,t) sc.E(r,t).*r, sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.ya=1/sc.EA*integral2(@(r,t) r.*cos(t).*r.*sc.E(r,t),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.zcg=1/sc.EA*integral2(@(r,t) r.*sin(t).*r.*sc.E(r,t),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.EJy=integral2(@(r,tz) (r.*sin(t)).^2.*sc.E(y,z),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.EJz=integral2(@(r,tz) (r.*cos(t)).^2.*sc.E(y,z),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.EJyz=integral2(@(r,tz) r.*sin(t).*r.*cos(t).*sc.E(y,z),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    
else
    error('Section is not defined properly')
end
