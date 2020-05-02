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
function M = el_mass_assembly(m, zM, yM, Ix, Iy, Iz, Iyz, L)

% DoF organisation [u1 v1 w1 th1 phi1 psi1 u2 v2 w2 th2 phi2 psi2]'


M = [[        (L*m)/3,               (L*m*yM)/4,                (L*m*zM)/4,              0,             (L^2*m*zM)/24,            -(L^2*m*yM)/24,        (L*m)/6,               -(L*m*yM)/4,              -(L*m*zM)/4,              0,            -(L^2*m*zM)/24,             (L^2*m*yM)/24];
[     (L*m*yM)/4,    (L*(21*Iz + 26*m))/70,              (3*Iyz*L)/10, -(7*L*m*zM)/20,             -(Iyz*L^2)/40,  (L^2*(21*Iz + 44*m))/840,     (L*m*yM)/4,    -(3*L*(7*Iz - 3*m))/70,            -(3*Iyz*L)/10, -(3*L*m*zM)/20,             -(Iyz*L^2)/40,  (L^2*(21*Iz - 26*m))/840];
[     (L*m*zM)/4,             (3*Iyz*L)/10,     (L*(21*Iy + 26*m))/70,  (7*L*m*yM)/20, -(L^2*(21*Iy + 44*m))/840,              (Iyz*L^2)/40,     (L*m*zM)/4,             -(3*Iyz*L)/10,   -(3*L*(7*Iy - 3*m))/70,  (3*L*m*yM)/20, -(L^2*(21*Iy - 26*m))/840,              (Iyz*L^2)/40];
[              0,           -(7*L*m*zM)/20,             (7*L*m*yM)/20,       (Ix*L)/3,            -(L^2*m*yM)/20,            -(L^2*m*zM)/20,              0,            -(3*L*m*zM)/20,            (3*L*m*yM)/20,       (Ix*L)/6,             (L^2*m*yM)/30,             (L^2*m*zM)/30];
[  (L^2*m*zM)/24,            -(Iyz*L^2)/40, -(L^2*(21*Iy + 44*m))/840, -(L^2*m*yM)/20,    (L^3*(7*Iy + 2*m))/210,             -(Iyz*L^3)/30, -(L^2*m*zM)/24,              (Iyz*L^2)/40, (L^2*(21*Iy - 26*m))/840, -(L^2*m*yM)/30,   -(L^3*(7*Iy + 6*m))/840,             (Iyz*L^3)/120];
[ -(L^2*m*yM)/24, (L^2*(21*Iz + 44*m))/840,              (Iyz*L^2)/40, -(L^2*m*zM)/20,             -(Iyz*L^3)/30,    (L^3*(7*Iz + 2*m))/210,  (L^2*m*yM)/24, -(L^2*(21*Iz - 26*m))/840,            -(Iyz*L^2)/40, -(L^2*m*zM)/30,             (Iyz*L^3)/120,   -(L^3*(7*Iz + 6*m))/840];
[        (L*m)/6,               (L*m*yM)/4,                (L*m*zM)/4,              0,            -(L^2*m*zM)/24,             (L^2*m*yM)/24,        (L*m)/3,               -(L*m*yM)/4,              -(L*m*zM)/4,              0,             (L^2*m*zM)/24,            -(L^2*m*yM)/24];
[    -(L*m*yM)/4,   -(3*L*(7*Iz - 3*m))/70,             -(3*Iyz*L)/10, -(3*L*m*zM)/20,              (Iyz*L^2)/40, -(L^2*(21*Iz - 26*m))/840,    -(L*m*yM)/4,     (L*(21*Iz + 26*m))/70,             (3*Iyz*L)/10, -(7*L*m*zM)/20,              (Iyz*L^2)/40, -(L^2*(21*Iz + 44*m))/840];
[    -(L*m*zM)/4,            -(3*Iyz*L)/10,    -(3*L*(7*Iy - 3*m))/70,  (3*L*m*yM)/20,  (L^2*(21*Iy - 26*m))/840,             -(Iyz*L^2)/40,    -(L*m*zM)/4,              (3*Iyz*L)/10,    (L*(21*Iy + 26*m))/70,  (7*L*m*yM)/20,  (L^2*(21*Iy + 44*m))/840,             -(Iyz*L^2)/40];
[              0,           -(3*L*m*zM)/20,             (3*L*m*yM)/20,       (Ix*L)/6,            -(L^2*m*yM)/30,            -(L^2*m*zM)/30,              0,            -(7*L*m*zM)/20,            (7*L*m*yM)/20,       (Ix*L)/3,             (L^2*m*yM)/20,             (L^2*m*zM)/20];
[ -(L^2*m*zM)/24,            -(Iyz*L^2)/40, -(L^2*(21*Iy - 26*m))/840,  (L^2*m*yM)/30,   -(L^3*(7*Iy + 6*m))/840,             (Iyz*L^3)/120,  (L^2*m*zM)/24,              (Iyz*L^2)/40, (L^2*(21*Iy + 44*m))/840,  (L^2*m*yM)/20,    (L^3*(7*Iy + 2*m))/210,             -(Iyz*L^3)/30];
[  (L^2*m*yM)/24, (L^2*(21*Iz - 26*m))/840,              (Iyz*L^2)/40,  (L^2*m*zM)/30,             (Iyz*L^3)/120,   -(L^3*(7*Iz + 6*m))/840, -(L^2*m*yM)/24, -(L^2*(21*Iz + 44*m))/840,            -(Iyz*L^2)/40,  (L^2*m*zM)/20,             -(Iyz*L^3)/30,    (L^3*(7*Iz + 2*m))/210]];
 
end