clear all, clc
syms lambda q c Vinf e L x real 

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 

Ca = zeros(12)*x; 

A1 = zeros(4)*x; 
A1(4,4) = (-q*c^2*pi/4*c/(2*Vinf))*cos(lambda)^2 ...
    + e*cos(lambda)^2*2*pi*c/(2*Vinf)*q*c; %...
    %-e*cos(lambda)*2*pi*e*cos(lambda)/Vinf;
Ca = Ca + int(N'*A1*N,x,-1,1)*L/2*cos(lambda); 

A2 = zeros(4)*x; 
A2(4,3) = q*c^2*pi/4*c/(2*Vinf)*cos(lambda)*sin(lambda) ...
    - e*cos(lambda)*2*pi*c/(2*Vinf)*sin(lambda)*q*c; %...
    %+ e*cos(lambda)*2*pi*e*sin(lambda)/Vinf; 
A2 = -A2; %sx / th
Ca = Ca + int(N'*A2*diff(N),x,-1,1)*cos(lambda);

A3 = zeros(4)*x; 
A3(3,4) = q*c^2*pi/4*c/(2*Vinf)*cos(lambda)*sin(lambda) ...
    -e*2*pi*c/(2*Vinf)*sin(lambda)*cos(lambda)*q*c; % ...
    %+ e*sin(lambda)*2*pi*e*cos(lambda)/Vinf; 
A3 = -A3; %sx / th
Ca = Ca + int(diff(N')*A3*N, x, -1,1)*cos(lambda); 

A4 = zeros(4)*x; 
A4(3,3) = -q*c^2*pi/4*c/(2*Vinf)*sin(lambda)^2 ...
    + e*sin(lambda)^2*2*pi*c/(2*Vinf)*q*c; %... 
    %-e*sin(lambda)*2*pi*e*sin(lambda)/Vinf; 
Ca = Ca + int(diff(N')*A4*diff(N),x,-1,1)*2/L*cos(lambda); 

A5 = zeros(4)*x; 
A5(3,3) = -2*pi/Vinf*q*c; 
Ca = Ca + int(N'*A5*N,x,-1,1)*L/2*cos(lambda);

A6 = zeros(4)*x; 
A6(3,4) = 2*pi*c/(2*Vinf)*cos(lambda)*q*c; %-2*pi*e*cos(lambda)/Vinf; 
A6 = -A6; %sx / th
Ca = Ca + int(N'*A6*N,x,-1,1)*L/2*cos(lambda); 

A7 = zeros(4)*x; 
A7(3,3) = -2*pi*c/(2*Vinf)*sin(lambda)*q*c; %+2*pi*e*sin(lambda)/Vinf;
Ca = Ca + int(N'*A7*diff(N),x,-1,1)*cos(lambda); 

A8 = zeros(4)*x; 
A8(4,3) = -e*cos(lambda)*2*pi/Vinf*q*c; 
A8 = -A8; %sx / th
Ca = Ca + int(N'*A8*N,x,-1,1)*L/2*cos(lambda); 

q = 1; 
Vinf = 1; 
eval(Ca)