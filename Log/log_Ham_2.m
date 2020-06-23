%% bla bla 
clear all, clc, close all
syms Q11 Q12 Q13 Q21 Q22 Q23 Q31 Q32 Q33 L x m k F G b e chord lambda
% m = el.sc.m
% b = chord/2
%% Define geometrical variables
a = -1/2; % we want to be in the aerodynamic centre
% c = 1/2;  % we put the aileron hinge at 75% of the chord

F = real(besselh(1, 2, k)./(besselh(1, 2, k) + 1i*besselh(0, 2, k)));
G = imag(besselh(1, 2, k)./(besselh(1, 2, k) + 1i*besselh(0, 2, k)));

% h, alpha 
coeff(1,1) = pi*(k^2-2*1i*k*(F+1i*G))*2/chord; 
coeff(1,2) = pi*(1i*k+k^2*a+2*(1+1i*k*(1/2-a))*(F+1i*G));
coeff(2,1) = -1i*k*pi/2*(F+1i*G)*2/chord; 
coeff(2,2) = pi/2*(-1i*k*1/2+k^2*1/8+(1+1i*k*(1/2-a))*(F+1i*G)); 

% want to pass into coordinates u,v,w,th
T = [0 0 1 e*cos(lambda) -e*sin(lambda); 
    0 0 0 cos(lambda) -sin(lambda)];

dim = [chord, 0; 0, chord^2];
coeff = transpose(T)*dim*coeff*T; 

% Original shape functions
N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0 ;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0];
% Added shape function for w/x
N = [N; diff(N(3,:),x)*2/L]; 

Ham = int(transpose(N)*coeff*N,x,-1,1)*L/2*cos(lambda); 

