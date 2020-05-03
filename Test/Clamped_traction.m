% This is a test:
% Displacement of a clamped beam subject to pure traction
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
N=1e3;          % Traction

% /|
% /|------------------------  ----> N
% /|
%% Build beam

beam=b_constant_p_square(L,l,E,G,rho,nel);

%% Constraints
K=beam.K(7:end,7:end);
M=beam.M(7:end,7:end);

%% Force

f=zeros(size(M,1),1);
f(end-5)=N;

%% Exact solution

exact_u_right=N/beam.el(1).sc.EA.*beam.L;

%% Solution

u=f\K;
u_right=u(end-5);

%% Disp
disp('Calculated u:')
disp(u_right)
disp('Exact u:')
disp(exact_u_right)
