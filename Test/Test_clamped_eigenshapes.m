% This is a test:
% Transverse vibrations of clamped-free (cantilever) uniform Euler-Bernoulli beam).
cd ..
init
clear all, close all, clc

%% Data
L=1e3;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=20;         % Number of elements
l=10;           % Side of the square
N=3;           % Number of eigenshapes i want to evaluate

%% Build beam

beam=b_constant_p_square_test(L,l,E,G,rho,nel);

%% Constraints
K=beam.K(7:end,7:end);
M=beam.M(7:end,7:end);

%% Solution
plot_cond=0;
w = ROM_solver(N, M, K, nel, L, plot_cond);

%% Exact solution

exact_k=[1.875104 4.694091 7.854757 10.995541 14.137168]';
exact_w=sqrt(beam.el(1).sc.EJy/beam.el(1).sc.m/L^4).*exact_k.^2;

%% Disp
disp(w)
disp(exact_w)


