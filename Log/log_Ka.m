clear all, clc, close all

syms x lambda c e CLa CMac CL0 L real

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 

Ka = zeros(12)*x; 

A1 = zeros(4)*x; 
A1(4,4) = cos(lambda)^3*e*c*CLa; 
Ka = Ka + int(N'*A1*N,x,-1,1)*L/2; 

A2 = zeros(4)*x; 
A2(4,3) = -cos(lambda)^2*sin(lambda)*e*c*CLa;
%A2 = -A2; %dx / w/x
A2 = -A2; %sx / th
Ka = Ka + int(N'*A2*diff(N,x),x,-1,1);

A3 = zeros(4)*x; 
A3(3,4) = cos(lambda)^2*c*CLa; 
A3 = -A3; %sx / th
A4 = zeros(4)*x; 
A4(3,4) = -cos(lambda)^2*sin(lambda)*e*c*CLa; 
%A4 = - A4; %dx w/x
A4 = -A4; %sx / th
Ka = Ka + int(N'*A3*N,x,-1,1)*L/2 + int(diff(N',x)*A4*N,x,-1,1);

A5 = zeros(4)*x; 
A5(3,3) = -cos(lambda)*sin(lambda)*c*CLa;
%A5 = - A5; %dx w/x
A6 = zeros(4)*x; 
A6(3,3) = cos(lambda)*sin(lambda)^2*e*c*CLa;
Ka = Ka + int(N'*A5*diff(N,x),x,-1,1) + int(diff(N',x)*A6*diff(N,x),x,-1,1)*2/L;

Ka(:,4+6) = -Ka(:,4+6); 
Ka(:,5+6) = -Ka(:,5+6); 
Ka(4+6,:) = -Ka(4+6,:); 
Ka(5+6,:) = -Ka(5+6,:); 
