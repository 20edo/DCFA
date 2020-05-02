function [sc]=sc_inertia(sc)

% Calculate the inertia propreties of a section
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
    sc.m=integral2(sc.rho,sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.ycg=1/sc.m*integral2(@(y,z) y.*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.zcg=1/sc.m*integral2(@(y,z) z.*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.Iy=integral2(@(y,z) z.^2.*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.Iz=integral2(@(y,z) y.^2.*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.Iyz=integral2(@(y,z) y.*z.*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.Jp=integral2(@(y,z) (y.^2+z.^2).*sc.rho(y,z),sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
%% Polar coordinates        y=r.*cos(t)/z=r.sin(t)
elseif sc.pol==true
    sc.m=integral2(@(t,r) sc.rho(r,t), sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.ycg=1/sc.m*integral2(@(t,r) r.*cos(t).*r.*sc.rho(r,t),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.zcg=1/sc.m*integral2(@(t,r) r.*sin(t).*r.*sc.rho(r,t),sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.Iy=integral2(@(t,r) (r.*sin(t)).^2.*sc.rho(r,t).*r,sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.Iz=integral2(@(t,r) (r.*cos(t)).^2.*sc.rho(r,t).*r,sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.Iyz=integral2(@(t,r) r.*sin(t).*r.*cos(t).*sc.rho(r,t).*r,sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.Jp=integral2(@(t,r) r.^2.*sc.rho(r,t).*r,sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
else
    error('Section is not defined properly')
end
