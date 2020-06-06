syms x lambda c e CLa CMac CL0 L real

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 

fa = zeros(12,1)*x; 

A1 = zeros(4,1)*x; 
A1(4,1) = cos(lambda)^2*(c^2*CMac + e*c*CL0); 
A1 = -A1; %sx
fa = fa + int(N'*A1,x,-1,1)*L/2; 

A2 = zeros(4,1)*x; 
A2(3,1) = cos(lambda)*c*CL0; 
fa = fa + int(N'*A2,x,-1,1)*L/2; 

A3 = zeros(4,1)*x; 
A3(3,1) = -cos(lambda)*sin(lambda)*(c^2*CMac + e*c*CL0); 
% A3 = -A3; %dx
fa = fa + int(diff(N')*A3,x,-1,1); 

