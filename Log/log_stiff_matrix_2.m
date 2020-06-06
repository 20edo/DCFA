clear all 
syms EA z_A y_A GJ EJy EJz EJyz s ell real
% moltiplicando le matrici che fanno cinematica
% le proprietà inerziali sono definite a partire dall'origine degli assi
syms z y u v w th real

N = [(1-s)/2   0   0   0   0   0   (1+s)/2   0   0   0   0   0;
    0, 1/4*(2-3*s+s^3), 0, 0, 0, 1/4*(1-s-s^2+s^3)*ell/2, 0, 1/4*(2+3*s-s^3), 0,0,0,1/4*(-1-s+s^2+s^3)*ell/2;
    0,0,1/4*(2-3*s+s^3),0,-1/4*(1-s-s^2+s^3)*ell/2,0,0,0,1/4*(2+3*s-s^3),0,-1/4*(-1-s+s^2+s^3)*ell/2,0; 
    0,0,0,(1-s)/2,0,0,0,0,0,(1+s)/2,0,0]; 
B = [1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1]; 
BN = B*N; 

BN(1,:) = diff(BN(1,:),s)*2/ell;
BN(2,:) = diff(diff(BN(2,:),s),s)*(2/ell)^2;
BN(3,:) = diff(diff(BN(3,:),s),s)*(2/ell)^2;
BN(4,:) = diff(BN(4,:),s)*2/ell;

% D = [EA, EA*zA -EA*yA 0; 
%     EA*zA EJy EJyz 0; 
%     -EA*yA EJyz EJz 0; 
%     0 0 0 GJ]; 


D = [EA, -EA*y_A -EA*z_A 0; 
    -EA*y_A EJz EJyz 0; 
    -EA*z_A EJyz EJy 0; 
    0 0 0 GJ]; 

% vwp = [u v w th]*D*[u v w th]'; 




Stiff = simplify(int(BN'*D*BN*L/2,s,-1,1));
% 
% zA = 0; 
% yA = 0; 
% EJyz = 0; 

% eval(Stiff)