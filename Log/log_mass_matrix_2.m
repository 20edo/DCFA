clear all 
syms m z_M y_M Ix Iy Iz Iyz x ell real
% moelltipellicando elle matrici che fanno cinematica
% elle proprietà inerziaelli sono definite a partire daell centro di massa
syms z y real

N = [(1-x)/2   0   0   0   0   0   (1+x)/2   0   0   0   0   0;
    0, 1/4*(2-3*x+x^3), 0, 0, 0, 1/4*(1-x-x^2+x^3)*ell/2, 0, 1/4*(2+3*x-x^3), 0,0,0,1/4*(-1-x+x^2+x^3)*ell/2;
    0,0,1/4*(2-3*x+x^3),0,-1/4*(1-x-x^2+x^3)*ell/2,0,0,0,1/4*(2+3*x-x^3),0,-1/4*(-1-x+x^2+x^3)*ell/2,0; 
    0,0,0,(1-x)/2,0,0,0,0,0,(1+x)/2,0,0]; 
D = [1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1;
    0 0 0 0; 
    0 0 0 0]; 
D2 = [1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1;
    0 0 -1 0; 
    0 1 0 0]; 
DN = D*N; 
DN2 = D2*N; 
DN2(5,:) = diff(DN2(5,:)); 
DN2(6,:) = diff(DN2(6,:)); 
for i = 1:12
DN(5,i) = -diff(N(3,i));
DN(6,i) = diff(N(2,i));
end


T = [1 0 0, 0 z_M -y_M; 
    0 1 0, -z_M 0 0; 
    0 0 1, y_M 0 0; 
    0 0 0, 1 0 0; 
    0 0 0, 0 1 0; 
    0 0 0, 0 0 1]; 

TDN = T*DN; 
TDN2 = T*DN2; 

MMM = [m 0 0, 0 0 0; 
    0 m 0, 0 0 0; 
    0 0 m, 0 0 0; 
    0 0 0, Ix 0 0; 
    0 0 0, 0 Iy -Iyz; 
    0 0 0, 0 -Iyz Iz]; 

Mass = int(TDN'*MMM*TDN*ell/2,-1,1);
Mass2 = int(TDN2'*MMM*TDN2*ell/2,-1,1);

z_M = 0; 
y_M = 0; 
%Iy = 0; 
%Iz = 0; 
Iyz = 0; 

decap = eval(Mass);

