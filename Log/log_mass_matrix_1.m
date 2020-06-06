clear all
syms m zM yM Ix Iy Iz Iyz x L real
% moltiplicando le matrici che fanno cinematica
% le quantità inerziali sono definite a partire dal punto (0,0)
syms z y real

N = [1 0 0 0 z -y; 
    0 1 0 -z 0 0; 
    0 0 1 y 0 0; 
    0 0 0 1 0 0]; 




M = [m  0   0   0    m*zM   -m*yM; 
     0  m   0  -m*zM    0       0; 
     0  0    m  m*yM    0       0; 
     0 -m*zM m*yM Ix    0       0; 
   m*zM  0   0     0   Iy    -Iyz; 
   -m*yM 0   0     0  -Iyz    Iz]; 




% u     (1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0       u1
% v     0   0   0   0   0   0   0   0   0   0   0   0       v1
% w     0   0   0   0   0   0   0   0   0   0   0   0       w1
% th    0   0   0   0   0   0   0   0   0   0   0   0       th1
% phi   0   0   0   0   0   0   0   0   0   0   0   0       phi1
% psi   0   0   0   0   0   0   0   0   0   0   0   0       psi1    
%                                                           u2
%                                                           v2
%                                                           w2
%                                                           th2
%                                                           phi2
%                                                           psi2


N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0; 
    0,0,-1/4*(-3+3*x^2)*2/L,0,1/4*(-1-2*x+3*x^2),0,0,0,-1/4*(+3-3*x^2)*2/L,0,1/4*(-1+2*x+3*x^2),0; 
    0,1/4*(-3+3*x^2)*2/L,0,0,0,1/4*(-1-2*x+3*x^2),0,1/4*(+3-3*x^2)*2/L,0,0,0,1/4*(-1+2*x+3*x^2)];

Mass = int(N'*M*N*L/2,-1,1);

% zM = 0; 
% yM = 0; 
% Iy = 0; 
% Iz = 0; 
% Iyz = 0; 

disacc = eval(Mass);



