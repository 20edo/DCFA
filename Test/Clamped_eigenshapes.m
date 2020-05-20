% This is a test:
% Transverse vibrations of clamped-free (cantilever) uniform Euler-Bernoulli beam).
cd ..
init
clear all, clc

%% Data
L=100;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=10;         % Number of elements
l=10;           % Side of the square

%% Build beam
% if ~exist('gcp')
%     parpool(4);
% end
% 
% tic
% spmd
%     beam=b_constant_p_square_test(L,l,E,G,rho,nel);
% end
% t_parallel=toc;
% 
% t_serial=tic
% beam=b_constant_p_square_test(L,l,E,G,rho,nel);
% t_serial=toc(t_serial);

beam=b_constant_p_square(L,l,E,G,rho,nel);
%% Exact solution

exact_k=[1.875104 4.694091 7.854757 10.995541 14.137168]';
exact_w=sqrt(beam.el(1).sc.EJy/beam.el(1).sc.m/L^4).*exact_k.^2;

%% Constraints
K=beam.K(7:end,7:end);
M=beam.M(7:end,7:end);
%% Solution
[w, V] = ROM_solver(14, M, K);
w = real(w); 
[VV,D,FLAG]=eigs(K,M,100,'smallestabs');
w = diag(D).^0.5;
V = VV; 


%% Disp
disp(w)
disp(exact_w)