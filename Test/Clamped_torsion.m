% This is a test:
% Displacement of a clamped beam subject to shear in the y direction
cd ..
init
clear all, close all, clc

%% Data
L=100;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=5;          % Number of elements
l=10;           % Side of the square
Mx=3e7;         % Torque

%% Build beam

beam=b_constant_p_square(L,l,E,G,rho,nel);
for i=1:length(beam.el)
    beam.el(i).sc.ycg=0;
    beam.el(i).sc.zcg=0;
    beam.el(i).sc.ya=0;
    beam.el(i).sc.za=0;
    beam.el(i).sc.yct=0;
    beam.el(i).sc.yct=0;
end
%% Constraints
K=beam.K(7:end,7:end);
M=beam.M(7:end,7:end);

%% Force

f=zeros(size(M,1),1);
f(end-2)=Mx;

%% Exact solution

exact_th_right=Mx/beam.el(1).sc.GJ.*beam.L;

%% Solution

u=K\f;
th_right=u(end-2);

%% Disp
disp('Calculated th:')
disp(th_right)
disp('Exact th:')
disp(exact_th_right)
beam.name='ClampedTorsionBeam'

for i = 2:nel+1
    beam.in(i).d = u(1+6*(i-2):6*(i-1),1);
end
% options.plot        = 1;
% options.plotColor   = 'green';
% options.saveSTL     = 0;
% options.ovs         = 4;
% [fig] = b_plot3d(beam, options);


close all

torsion3d=m_init;
torsion3d.en=en_ground([0 0 0]');
beam.o=[0 0 0]';
beam.vx=[0 1 0]';
beam.vy=[1 0 0]';

torsion3d=m_add_beam(torsion3d,beam)
% 3dtorsion=m_compute_matrices(3dtorsion)


   
options.plot_original          = 1;
options.plot_deformed          = 1;
options.plotColor              = 'green';
options.saveSTL                = 0;
options.point_section          = 8;
[fig] = m_plot3d(torsion3d,options)

% return where I come from ... home sweet home
cd ./test
