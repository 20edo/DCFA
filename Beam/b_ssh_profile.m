function b=b_ssh_profile(L,c,h,t,E,G,rho,nel)
% Creates a beam with constant properties and tube section of side external radius R(x) and thickness t(x). 
% The length of the beam is L. nin is the number of internal nodes.
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
% Beam
b=b_init();
b.ssh=true;              % True if the section is defined as a profile (semishield)
b.c=@(x) c(x);              % Chord of the profile as a function of x
b.h=@(x) h(x);              % Heigth percentage of the profile as a function of x
b.t=@(x) t(x);              % Thickness of the profile as a funcion of x 

b.cart=true;
NACA = [0 0 1 2];
XX = str2num([num2str(NACA(3)),num2str(NACA(4))])/100;
b.Zmax=@(x,y) 5*XX.*c(x).*(0.2969*sqrt(y./c(x))-0.126*(y./c(x))-0.3516.*(y./c(x)).^2+0.2843.*(y./c(x)).^3-0.1015.*(y./c(x)).^4); 
b.Zmin=@(x,y) -5*XX.*c(x).*(0.2969*sqrt(y./c(x))-0.126*(y./c(x))-0.3516.*(y./c(x)).^2+0.2843.*(y./c(x)).^3-0.1015.*(y./c(x)).^4);
b.Ymax=@(x) c(x).*(1/2);
b.Ymin=@(x) -c(x).*(1/2);
% All functions must be defined in the same reference of the geometry (pol
% or cart)
b.E=@(x,y,z) E+0.*x+0.*y+0.*z;        % Young modulus of the point in the beam
b.G=@(x,y,z) G+0.*x+0.*y+0.*z;          % Shear modulus of the pooint in the beam
b.rho=@(x,y,z) rho+0.*x+0.*y+0.*z;       % Density of the point of the beam
% (Properties)
b.L=L;                   % Length of the beam     [m]
b.nel=nel;               % Number of elements

%% Speed up since the section is constant
b=b_build_elements_ssh_profile(b);
b=b_build_matrices(b);
end