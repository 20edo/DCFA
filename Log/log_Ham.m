%% bla bla 
clear all, clc, close all
syms Q11 Q12 Q13 Q21 Q22 Q23 Q31 Q32 Q33 L x m k F G b lambda
% m = el.sc.m
% b = chord/2
%% Define geometrical variables
a = -1/2; % we want to be in the aerodynamic centre
c = 1/2;  % we put the aileron hinge at 75% of the chord 
%% Define Theodorsen coefficients T
T1 = -1/3*sqrt(1-c^2)*(2+c^2)+c*acos(c); 
T2 = c*(1-c^2)-sqrt(1-c^2)*(1+c^2)*acos(c)+c*(acos(c))^2; 
T3 = -(1/8+c^2)*acos(c)^2+1/4*c*sqrt(1-c^2)*acos(c)*(7+2*c^2)-1/8*(1-c^2)*(5*c^2+4); 
T4 = -acos(c)+c*sqrt(1-c^2); 
T5 = -(1-c^2)-acos(c)^2+2*c*sqrt(1-c^2)*acos(c); 
T6 = T2; 
T7 = -(1/8+c^2)*acos(c)+1/8*c*sqrt(1-c^2)*(7+2*c^2); 
T8 = -1/3*sqrt(1-c^2)*(2*c^2+1)+c*acos(c); 
p = -1/3*(sqrt(1-c^2))^3; 
T9 = 1/2*(-p+a*T4); 
T10 = sqrt(1-c^2) + acos(c); 
T11 = acos(c)*(1-2*c)+sqrt(1-c^2)*(2-c); 
T12 = sqrt(1-c^2)*(2+c)-acos(c)*(2*c+1); 
T13 = 1/2*(-T7-(c-a)*T1); 
T14 = 1/16*1/2*a*c; 
%% Define Theodorsen functions
Q11r = 2*pi/m*(-(1/8+a^2)*k^2-2*(a^2-1/4)*G*k - 2*(a+1/2)*F);
Q11i = 2*pi/m*((1/2-a)*k + 2*(a^2-1/4)*F*k - 2*(a+1/2)*G); 
Q11 = Q11r + 1i*Q11i; 

Q21r = 2/m*((T7+(c-a)*T1)*k^2 -T12*(1/2-a)*G*k + T12*F); 
Q21i = 2/m*((p-T1-T4/2)*k+T12*(1/2-a)*F*k + T12*G); 
Q21 = Q21r + 1i*Q21i; 

Q31r = 2*pi/m*(a*k^2 - 2*(1/2-a)*G*k + 2*F); 
Q31i = 2*pi/m*(k+2*(1/2-a)*F*k + 2*G); 
Q31 = Q31r + 1i*Q31i; 

Q12r = 2/m*((T7+(c-a)*T1)*k^2+(T4+T10)+(a+1/2)*T11*G*k - 2*(a+1/2)*T10*F); 
Q12i = 2/m*(-(2*p+(1/2-a)*T4)*k-(a+1/2)*T11*F*k - 2(a+1/2)*T10*G); 
Q12 = Q12r + 1i*Q12i; 

Q22r = 2/(m*pi)*(T3*k^2 + (T5-T4*T10)-T11*T12/2*G*k*T10*T12*F); 
Q22i = 2/(m*pi)*(-T4*T11/2*k + T11*T12/2*F*k + T10*T12*G); 
Q22 = Q22r + 1i*Q22i; 

Q32r = 2/m*(T1*k^2 - T11*G*k + 2*T10*F); 
Q32i = 2/m*(-T4*k + T11*F*k + 2*T10*G); 
Q32 = Q32r + 1i*Q32i; 

Q13r = 2*pi/(M*b)*(a*k^2+2*(a+1/2)*G*k); 
Q13i = 2*pi/(M*b)*(-2*(a+1/2)*F*k); 
Q13 = Q13r + 1i*Q13i; 

Q23r = 2/(M*b)*(T1*k^2 - T12*G*k); 
Q23i = 2/(M*b)*(T12*F*k); 
Q23 = Q23r + 1i*Q23i; 

Q33r = 2*pi/(M*b)*(-k^2-2*G*k); 
Q33i = 2*pi/(M*b)*(2*F*k); 
Q33 = Q33r + 1i*Q33i; 


%% Transform the Q matrix from theodorsen reference to our reference
Q_th = [Q11, Q12, Q13;  % alpha, beta, h positive down
    Q21, Q22, Q23; 
    Q31, Q32, Q33]; 

% alpha = cos(lambda)*th - sin(lambda)*w/y 
% we want u,v,w,th,w/x,b positive up
T = [0 0 0 cos(lambda) sin(lambda) 0; 
    0 0 0 0 0 1; 
    0 0 -1 0 0 0]; 

Q_us = transpose(T)*Q_th*T; 

%% Expand the shape funcions

% Original shape functions
N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0 ;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0];
% Added shape function for w/x
N = [N; diff(N(3,:),x)*2/L]; 

N=  [N       zeros(size(N,1),1); 
    zeros(1,size(N,2)) 1];
%% Expand the Q in our reference through shape functions
Q_tot = transpose(N)*Q_us*N; 
Q_tot = simplify(Q_tot); 
Q_tot = int(Q_tot,x,-1,1)*L/2; 

Ham = Q_tot(1:end-1,1:end-1); 
Hbam = Q_tot(end,1:end-1); 
Hamb = Q_tot(1:end-1,end); 
Hbb = Q_tot(end,end); 










