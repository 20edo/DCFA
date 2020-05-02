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
function K = stiff_assembly_el(EA,zA,yA,GJ,EJy,EJz,EJyz,L)

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
end