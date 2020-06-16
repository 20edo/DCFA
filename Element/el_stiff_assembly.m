% This function gives the consisten element stiffness matrix given the
% elastic properties of the element (this is just the final assembly, 
% the reasoning behind is hidden)
%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Pasturenzi Lorenzo    944610
%               Tacchi Alberto        944579
%               Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
function el = el_stiff_assembly(el)

EA=el.sc.EA;
zA=el.sc.za;
yA=el.sc.ya;
GJ=el.sc.GJ;
EJy=el.sc.EJy;
EJz=el.sc.EJz;
EJyz=el.sc.EJyz;
L=el.L;

% DoF organisation [u1 v1 w1 th1 phi1 psi1 u2 v2 w2 th2 phi2 psi2]'

K = [[       EA/L,              0,              0,     0,     (EA*zA)/L,    -(EA*yA)/L,      -EA/L,              0,              0,     0,    -(EA*zA)/L,     (EA*yA)/L];
[          0,   (12*EJz)/L^3,  (12*EJyz)/L^3,     0, -(6*EJyz)/L^2,   (6*EJz)/L^2,          0,  -(12*EJz)/L^3, -(12*EJyz)/L^3,     0, -(6*EJyz)/L^2,   (6*EJz)/L^2];
[          0,  (12*EJyz)/L^3,   (12*EJy)/L^3,     0,  -(6*EJy)/L^2,  (6*EJyz)/L^2,          0, -(12*EJyz)/L^3,  -(12*EJy)/L^3,     0,  -(6*EJy)/L^2,  (6*EJyz)/L^2];
[          0,              0,              0,  GJ/L,             0,             0,          0,              0,              0, -GJ/L,             0,             0];
[  (EA*zA)/L,  -(6*EJyz)/L^2,   -(6*EJy)/L^2,     0,     (4*EJy)/L,   -(4*EJyz)/L, -(EA*zA)/L,   (6*EJyz)/L^2,    (6*EJy)/L^2,     0,     (2*EJy)/L,   -(2*EJyz)/L];
[ -(EA*yA)/L,    (6*EJz)/L^2,   (6*EJyz)/L^2,     0,   -(4*EJyz)/L,     (4*EJz)/L,  (EA*yA)/L,   -(6*EJz)/L^2,  -(6*EJyz)/L^2,     0,   -(2*EJyz)/L,     (2*EJz)/L];
[      -EA/L,              0,              0,     0,    -(EA*zA)/L,     (EA*yA)/L,       EA/L,              0,              0,     0,     (EA*zA)/L,    -(EA*yA)/L];
[          0,  -(12*EJz)/L^3, -(12*EJyz)/L^3,     0,  (6*EJyz)/L^2,  -(6*EJz)/L^2,          0,   (12*EJz)/L^3,  (12*EJyz)/L^3,     0,  (6*EJyz)/L^2,  -(6*EJz)/L^2];
[          0, -(12*EJyz)/L^3,  -(12*EJy)/L^3,     0,   (6*EJy)/L^2, -(6*EJyz)/L^2,          0,  (12*EJyz)/L^3,   (12*EJy)/L^3,     0,   (6*EJy)/L^2, -(6*EJyz)/L^2];
[          0,              0,              0, -GJ/L,             0,             0,          0,              0,              0,  GJ/L,             0,             0];
[ -(EA*zA)/L,  -(6*EJyz)/L^2,   -(6*EJy)/L^2,     0,     (2*EJy)/L,   -(2*EJyz)/L,  (EA*zA)/L,   (6*EJyz)/L^2,    (6*EJy)/L^2,     0,     (4*EJy)/L,   -(4*EJyz)/L];
[  (EA*yA)/L,    (6*EJz)/L^2,   (6*EJyz)/L^2,     0,   -(2*EJyz)/L,     (2*EJz)/L, -(EA*yA)/L,   -(6*EJz)/L^2,  -(6*EJyz)/L^2,     0,   -(4*EJyz)/L,     (4*EJz)/L]];



% Change of reference
R=[ 1   0   0   0   0   0
    0   1   0   -el.sc.zct  0  0
    0   0   1   el.sc.yct   0  0
    0   0   0   1   0   0
    0   0   0   0   1   0
    0   0   0   0   0   1];

R=sparse(R);

R=blkdiag(R,R);

K=R*K*transpose(R);


el.K=sparse(K);


el.K=K;
end