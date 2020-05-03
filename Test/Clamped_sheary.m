% This is a test:
% Displacement of a clamped beam subject to shear in the y direction
cd ..
init
clear all, close all, clc

%% Data
L=1e3;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=5;          % Number of elements
l=10;           % Side of the square
Ty=1e3;          % Traction

%% Build beam

beam=b_constant_p_square(L,l,E,G,rho,nel);

%% Constraints
K=beam.K(7:end,7:end);
M=beam.M(7:end,7:end);

%% Force

f=zeros(size(M,1),1);
f(end-3)=Ty;

%% Exact solution

exact_v_right=(-beam.L^3/6+beam.L^3/2)*Ty/beam.el(1).sc.EJy;

%% Solution

u=K\f;
v_right=u(end-3);

%% Disp
disp('Calculated u:')
disp(v_right)
disp('Exact u:')
disp(exact_v_right)
