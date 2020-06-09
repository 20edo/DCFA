%%

clear all, clc
syms lambda q c Vinf e L x a real 


Q = [q*c^2*(-pi/4)*c/(2*Vinf),                  0;    % alpha w
    q*c*2*pi*(1/2-a)*c/(2*Vinf),   -q*c*2*pi*(1/Vinf)];


% alpha = cos(lambda)*th - sin(lambda)*w/y 
% we want u,v,w,th,w/x positive up
% Move reference point to the elastic axis
% - h ca = -w -e*cos(lambda)*th + e*sin(lambda)*w/y
TT = [0 0 0 cos(lambda) -sin(lambda) ; 
    0 0 1 +e*cos(lambda) -e*sin(lambda) ]';
T = [0 0 0 cos(lambda) -sin(lambda) ; 
    0 0 1 0 0 ];

Q_us = TT*(Q)*T; 

%% Expand the shape funcions

% Original shape functions
N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0 ;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0];
% Added shape function for w/x
N = [N; diff(N(3,:),x)*2/L]; 

%% Expand the Q in our reference through shape functions
Ca = transpose(N)*Q_us*N; 
Ca = simplify(Ca); 
Ca = int(Ca,x,-1,1)*L/2*cos(lambda); 

a = -1/2;
q = 1; 
Vinf = 1; 
eval(Ca)

