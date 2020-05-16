% This function gives the consisten element mass matrix given the inertial
% properties of the element (this is just the final assembly, the reasoning
% behind is hidden)
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
function el = el_mass_assembly(el)

m=el.sc.m;
zM=el.sc.zcg;
yM=el.sc.ycg;
Ix=el.sc.Jp;
Iy=el.sc.Iy;
Iz=el.sc.Iz;
Iyz=el.sc.Iyz;
L=el.L;

% DoF organisation [u1 v1 w1 th1 phi1 psi1 u2 v2 w2 th2 phi2 psi2]'


M = [[ (L*m)/3,                   (m*yM)/2,                   (m*zM)/2,              0,               (L*m*zM)/12,              -(L*m*yM)/12,      (L*m)/6,                  -(m*yM)/2,                  -(m*zM)/2,              0,              -(L*m*zM)/12,               (L*m*yM)/12];
[     (m*yM)/2, (13*L*m)/35 + (6*Iz)/(5*L),              (6*Iyz)/(5*L), -(7*L*m*zM)/20,                   -Iyz/10,    (11*m*L^2)/210 + Iz/10,     (m*yM)/2,  (9*L*m)/70 - (6*Iz)/(5*L),             -(6*Iyz)/(5*L), -(3*L*m*zM)/20,                   -Iyz/10,    Iz/10 - (13*L^2*m)/420];
[     (m*zM)/2,              (6*Iyz)/(5*L), (13*L*m)/35 + (6*Iy)/(5*L),  (7*L*m*yM)/20,  - Iy/10 - (11*L^2*m)/210,                    Iyz/10,     (m*zM)/2,             -(6*Iyz)/(5*L),  (9*L*m)/70 - (6*Iy)/(5*L),  (3*L*m*yM)/20,    (13*m*L^2)/420 - Iy/10,                    Iyz/10];
[            0,             -(7*L*m*zM)/20,              (7*L*m*yM)/20,       (Ix*L)/3,            -(L^2*m*yM)/20,            -(L^2*m*zM)/20,            0,             -(3*L*m*zM)/20,              (3*L*m*yM)/20,       (Ix*L)/6,             (L^2*m*yM)/30,             (L^2*m*zM)/30];
[  (L*m*zM)/12,                    -Iyz/10,   - Iy/10 - (11*L^2*m)/210, -(L^2*m*yM)/20, (m*L^3)/105 + (2*Iy*L)/15,             -(2*Iyz*L)/15, -(L*m*zM)/12,                     Iyz/10,     Iy/10 - (13*L^2*m)/420, -(L^2*m*yM)/30, - (L^3*m)/140 - (Iy*L)/30,                (Iyz*L)/30];
[ -(L*m*yM)/12,     (11*m*L^2)/210 + Iz/10,                     Iyz/10, -(L^2*m*zM)/20,             -(2*Iyz*L)/15, (m*L^3)/105 + (2*Iz*L)/15,  (L*m*yM)/12,     (13*m*L^2)/420 - Iz/10,                    -Iyz/10, -(L^2*m*zM)/30,                (Iyz*L)/30, - (L^3*m)/140 - (Iz*L)/30];
[      (L*m)/6,                   (m*yM)/2,                   (m*zM)/2,              0,              -(L*m*zM)/12,               (L*m*yM)/12,      (L*m)/3,                  -(m*yM)/2,                  -(m*zM)/2,              0,               (L*m*zM)/12,              -(L*m*yM)/12];
[    -(m*yM)/2,  (9*L*m)/70 - (6*Iz)/(5*L),             -(6*Iyz)/(5*L), -(3*L*m*zM)/20,                    Iyz/10,    (13*m*L^2)/420 - Iz/10,    -(m*yM)/2, (13*L*m)/35 + (6*Iz)/(5*L),              (6*Iyz)/(5*L), -(7*L*m*zM)/20,                    Iyz/10,  - Iz/10 - (11*L^2*m)/210];
[    -(m*zM)/2,             -(6*Iyz)/(5*L),  (9*L*m)/70 - (6*Iy)/(5*L),  (3*L*m*yM)/20,    Iy/10 - (13*L^2*m)/420,                   -Iyz/10,    -(m*zM)/2,              (6*Iyz)/(5*L), (13*L*m)/35 + (6*Iy)/(5*L),  (7*L*m*yM)/20,    (11*m*L^2)/210 + Iy/10,                   -Iyz/10];
[            0,             -(3*L*m*zM)/20,              (3*L*m*yM)/20,       (Ix*L)/6,            -(L^2*m*yM)/30,            -(L^2*m*zM)/30,            0,             -(7*L*m*zM)/20,              (7*L*m*yM)/20,       (Ix*L)/3,             (L^2*m*yM)/20,             (L^2*m*zM)/20];
[ -(L*m*zM)/12,                    -Iyz/10,     (13*m*L^2)/420 - Iy/10,  (L^2*m*yM)/30, - (L^3*m)/140 - (Iy*L)/30,                (Iyz*L)/30,  (L*m*zM)/12,                     Iyz/10,     (11*m*L^2)/210 + Iy/10,  (L^2*m*yM)/20, (m*L^3)/105 + (2*Iy*L)/15,             -(2*Iyz*L)/15];
[  (L*m*yM)/12,     Iz/10 - (13*L^2*m)/420,                     Iyz/10,  (L^2*m*zM)/30,                (Iyz*L)/30, - (L^3*m)/140 - (Iz*L)/30, -(L*m*yM)/12,   - Iz/10 - (11*L^2*m)/210,                    -Iyz/10,  (L^2*m*zM)/20,             -(2*Iyz*L)/15, (m*L^3)/105 + (2*Iz*L)/15]];
 

% Change of reference

R=[ 1   0   0   0   0   0
    0   1   0   el.sc.yct   0  0
    0   0   1   -el.sc.zct  0  0
    0   0   0   1   0   0
    0   0   0   0   1   0
    0   0   0   0   0   1];

R=sparse(R);

R=blkdiag(R,R);

M=R*M*transpose(R);


el.M=sparse(M);

el.M=M;
end