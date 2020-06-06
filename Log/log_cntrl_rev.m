clear all
syms x lambda c e CLa L CMb CLb real

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*L/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*L/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*L/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*L/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 

fb = zeros(12,1)*x; 

A1 = zeros(4,1)*x; 
A1(4,1) = cos(lambda)^2*(c^2*CMb + e*c*CLb); 
%A1 = -A1; %sx / th 
fb = fb + int(N'*A1,x,-1,1)*L/2; 

A2 = zeros(4,1)*x; 
A2(3,1) = cos(lambda)*c*CLb; 
fb = fb + int(N'*A2,x,-1,1)*L/2; 

A3 = zeros(4,1)*x; 
A3(3,1) = -cos(lambda)*sin(lambda)*(c^2*CMb + e*c*CLb); 
%A3 = -A3; %dx / w/x 
fb = fb + int(diff(N')*A3,x,-1,1); 

Lq = zeros(1,12)*x; 

A1 = zeros(1,4)*x; 
A1(1,3) = -c*CLa*sin(lambda)*cos(lambda);
%A1 = -A1; %dx /w/x 
Lq = Lq + int(A1*diff(N),x,-1,1); 

A2 = zeros(1,4)*x; 
A2(1,4) = c*CLa*cos(lambda)^2; 
%A2 = - A2; %sx / th 
Lq = Lq + int(A2*N,x,-1,1)*L/2; 
 
Lb = int(c*CLb*cos(lambda),x,-1,1)*L/2;



