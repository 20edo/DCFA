function sc=sc_ssh_profile(c,h,A,t,E,G,rho)

% Returns a profile section. The geometry is based on a NACA00XX whereas
% the properties are based on semi-shield approximation.
% c         Chord of the profile
% h         Heigth percentage of the profile (ex. NACA 0012->h 0.12)
% A         Total area of the section (panels+correnti)
% t         Thickness of the panels

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

NACA = [0 0 1 2];

% section 
sc=sc_init();

sc.pol=false;               % True if the section is defined in polar coordinates
sc.Rhomin=@(th) (R-t)+0.*th;        % Coordinate of the 'lower' rho boundary
sc.Rhomax=@(th) R+0.*th;        % Coordinate of the 'upper' rho boundary
sc.Thmin=0;               % Coordinate of the 'lower' th boundary
sc.Thmax=2*pi;               % Coordinate of the 'upper' th boundary
sc.E=@(Rho,th) E+0.*Rho+0.*th;
sc.G=@(Rho,th) G+0.*Rho+0.*th;
sc.rho=@(Rho,th) rho+0.*Rho+0.*th;

sc.cart=true;
XX = str2num([num2str(NACA(3)),num2str(NACA(4))])/100;
sc.Zmax=@(y) 5*XX*c*(0.2969*sqrt(y/c)-0.126*(y/c)-0.3516*(y/c).^2+0.2843*(y/c).^3-0.1015*(y/c).^4); 
sc.Zmin=@(y) -5*XX*c*(0.2969*sqrt(y/c)-0.126*(y/c)-0.3516*(y/c).^2+0.2843*(y/c).^3-0.1015*(y/c).^4);
sc.Ymax= c/2;
sc.Ymin= -c/2;
sc.E=@(y,z) E+0.*y+0.*z;
sc.G=@(y,z) G+0.*y+0.*z;
sc.rho=@(y,z) rho+0.*y+0.*z;
% sc=sc_compute_property(sc);


A=A/8;

%% Solution of the section (xy coordinates)
% (To understand see hand calculations)

% Shear center problem

sc.zct=0;   % Simmetry wrt y axis 
sc.yct=0;   % Imposed

h=c*h;
xcg=5/16*c;
L=sqrt((1/4)^2+(h/2)^2);

EQUATIONS= [    2*L+h/3     -h/3    0           % th1=0
                -h/3    c+2*L+h/3   0           % th2=0
                1           5      4/c/h ];     % Torque equilibrium CG
    
b=[ 0
    2*c/3/h+2*L/h
    49/12/h];

% x=[qast1 quast2 d]' Unknown vector
x=EQUATIONS\b;
d=x(3);

xct= xcg+d;

% Shear stifness 

EQUATIONS=[ 10*L+2*h            -c-2*h-2*L          0            % th1=th2
            xct*h/2+c*h/8       9/8*c*h-h/2*xct     0            % Torque equilibrium CG
            2*L+h/3             -h/3                -G*c*h*t/4]; % th calculation

b=[ 0
    1/2
    0];

% x= [quast1 quast2 th] unknown vector
x=EQUATIONS\b;
th=x(3);

sc.GJ=1/th;

%% Mass properties

sc.ycg=-d;
sc.zcg=0;       % Simmetry
sc.EA=8*A*E;
sc.GA=8*A*G;
sc.Iy=3/2*A*h^2*rho;
sc.Iz=xct^2*A+(c/4-xct)^2*4*A+(3/4*c-xct)^2*2*A+(c-xct)^2*A; 
sc.Iyz=0;       % Simmetry
sc.Jp=xct^2*A+((c/4-xct)^2+(h/2)^2)*4*A+((c*3/4-xct)^2+(h/2)^2)*2*A+(c-xct)^2*A;
sc.m=8*A*rho;


%% Elastic properties

sc.ya=-d;
sc.za=0;
sc.EJy=E*sc.Iy;
sc.EJz=E*sc.Iz;
sc.EJyz=0;
%% Keep track of the origin of the geometry
sc.yo=-xct;
end
