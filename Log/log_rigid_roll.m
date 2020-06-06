%% Bla Bla 
clear all; clc; close all; 
syms Jx Ja m lambda x x1 x2 c CLa CLb L d y e real

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 

y = ([(1-x)/2 (1+x)/2]*[x1,x2]')*cos(lambda); 

Jx = 2*int(m*y^2*cos(lambda),x,-1,1)*L/2;
lp = -2*int(y^2*c*CLa*cos(lambda),x,-1,1)*L/2;
lb = 2*int(y*c*CLb,x,-1,1)*L/2;

lq = zeros(1,12)*x; 

A1 = zeros(1,4)*x; 
A1(1,4) = 2*(y*c*CLa*cos(lambda)*cos(lambda)); 
lq = lq + int(A1*N,x,-1,1)*L/2; 

A2 = zeros(1,4)*x; 
A2(1,3) = -2*(y*c*CLa*sin(lambda)*cos(lambda)); 
lq = lq + int(A2*diff(N),x,-1,1); 

clear A1
clear A2 

Sq = zeros(12,1)*x; 

A1 = zeros(4,1)*x; 
A1(4,1) = -m*y*d*cos(lambda); 
Sq = Sq + int(N'*A1,x,-1,1)*L/2; 

A2 = zeros(4,1)*x; 
A2(3,1) = +m*y*cos(lambda);
Sq = Sq + int(N'*A2,x,-1,1)*L/2; 

clear A1 
clear A2 

fp = zeros(12,1)*x; 

A1 = zeros(4,1)*x; 
A1(4,1) = -y*e*c*CLa*cos(lambda); 
fp = fp + int(N'*A1,x,-1,1)*L/2; 

A2 = zeros(4,1)*x; 
A2(3,1) = -y*c*CLa*cos(lambda);
fp = fp + int(N'*A2,x,-1,1)*L/2; 










