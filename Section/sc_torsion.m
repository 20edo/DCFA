function sc=sc_torsion(sc,Hmax,np)
% Calculate torsion factor of the section by a FEM model, with a
% displacement approach.
% Check for sc.pol==true
% sc        ->      Section
% Hmax      ->      Max characteristic dimension of the FEM triangle    [optional]
% np        ->      Number of points used to define the geometry        [optional]
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

%% Define missing input

if nargin==1
    np=1e2;
    if sc.cart==true
        Hmax=sqrt(sc.Ymax^2+sc.Zmax(sc.Ymax)^2)/10;
    else
        Hmax=sc.Rhomax(sc.Thmin)/10;
    end
elseif nargin ==2
    np=1e2;
end
    
%% Define geometry

if sc.cart==true
    % Discretize boundary
    y=[linspace(sc.Ymin, sc.Ymax, np) linspace(sc.Ymax, sc.Ymin, np)];
    z=[sc.Zmin(linspace(sc.Ymin, sc.Ymax, np)) sc.Zmax(linspace(sc.Ymax, sc.Ymin, np))];
    % Eliminate duplicate points (Check, not working)
    %y=y([2:end/2, end/2+2:end]);
    %z=z([2:end/2, end/2+2:end]);
    pgon=polyshape(y,z);
    Edges=[1:4];
elseif sc.pol==true
    th=linspace(sc.Thmin, sc.Thmax, np);
    pgon=polyshape({sc.Rhomax(th).*cos(th),sc.Rhomin(th).*cos(th)},...
        {sc.Rhomax(th).*sin(th),sc.Rhomin(th).*sin(th)});
    Edges=[1];
    y=0;
    z=0;
else
    error('Section is not defined properly')
end
    
model = createpde;

%% Mesh

tr=triangulation(pgon);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
%pdegplot(model,'EdgeLabels','on')
generateMesh(model,'Hmax',Hmax);
%pdemesh(model) % Plot della mesh
%% Scrivo l'equazione (Laplaciano)
specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',0); 
%% B.C
bc=@(location,state) location.y.*location.nx-location.x.*location.ny;
applyBoundaryCondition(model,'neumann','Edge',Edges,'g',bc);

%% Solution
result = solvepde(model);

if sc.cart==true
    sc.GJ=integral2(@(y,z) integrand(y,z,sc,result),...
        sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax);
    sc.yct=-1/sc.Iy*integral2(@(y,z) integrand2(y,z,sc,result),...
        sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax,...
        'AbsTol',1e-3);
    sc.zct=1/sc.Iy*integral2(@(y,z) integrand3(y,z,sc,result),...
    sc.Ymin,sc.Ymax,sc.Zmin,sc.Zmax,...
        'AbsTol',1e-3);
elseif sc.pol==true
    sc.GJ=integral2(@(t,r) integrand(y,z,sc,result,r,t),...
        sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.yct=-1/sc.Iy*integral2(@(t,r) integrand2(y,z,sc,result,r,t),...
        sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
    sc.zct=1/sc.Iz*integral2(@(t,r) integrand3(y,z,sc,result,r,t),...
        sc.Thmin,sc.Thmax,sc.Rhomin,sc.Rhomax);
end
