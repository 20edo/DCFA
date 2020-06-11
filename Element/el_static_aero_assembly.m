% Function that gives the expression of the static aerodynamic matrices
%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%
%
%

function el=el_static_aero_assembly(b,el,lambda,i,alpha,dx)

%% Assign geometrical variables

L=el.L;                         % Length of the element
c = el.sc.Ymax - el.sc.Ymin ;   % Chord
e = -el.sc.yo - 1/4*c;          % e -> distance between the shear center (elastic axis) and AC
d = el.sc.ycg;                  % distance of the mass center wrt shear center
m = el.sc.m;                    % mass of the element

pos1 = b.o + b.vx*b.in(i).x;    % Position of the first node of the element
pos2 = b.o + b.vx*b.in(i+1).x;  % Position of the second node of the element
x1 = pos1(2);
x2 = pos2(2);

%% Assign the aerodinamic coefficients

% Profile coefficients
CLa = el.sc.CLa;                % CL/alpha
CL0 = alpha*CLa;                % CL at zero incidence
CMac = el.sc.CMaca*(-alpha);    % Moment coefficient wrt aerodynamic center

% Flap coefficients
CLb=el.sc.CLb;                  % CL/beta
CMb=el.sc.CMb;                  % CM/beta wrrt aerodynamic center

%% Right wing

if dx==1
    % Aerodynamic stiffness matrix
    Ka = [[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]
[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]
[ 0, 0, (CLa*c*sin(2*lambda))/4 - (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),    (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2 - (7*CLa*L*c*cos(lambda)^2)/20,      (CLa*L*c*sin(2*lambda))/20 + (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0, 0, 0,   (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L) - (CLa*c*sin(2*lambda))/4,    (3*CLa*L*c*cos(lambda)^2)/20 - (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,        (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0]
[ 0, 0,                                 (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                   (CLa*L*c*e*cos(lambda)^3)/3,                                 (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                                  -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                  -(CLa*L*c*e*cos(lambda)^3)/6,                                   (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
[ 0, 0,   (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10 - (CLa*L*c*sin(2*lambda))/20, (CLa*L^2*c*cos(lambda)^2)/20 + (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12,                              -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0, 0, 0,     (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12 - (CLa*L^2*c*cos(lambda)^2)/30, - (CLa*L^2*c*sin(2*lambda))/120 - (CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/30, 0]
[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]
[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]
[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]
[ 0, 0, (CLa*c*sin(2*lambda))/4 + (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),  - (3*CLa*L*c*cos(lambda)^2)/20 - (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,    - (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0, 0, 0, - (CLa*c*sin(2*lambda))/4 - (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),    (7*CLa*L*c*cos(lambda)^2)/20 + (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,        (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10 - (CLa*L*c*sin(2*lambda))/20, 0]
[ 0, 0,                                -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                  -(CLa*L*c*e*cos(lambda)^3)/6,                                 (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                                   (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                   (CLa*L*c*e*cos(lambda)^3)/3,                                   (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
[ 0, 0, - (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, (CLa*L^2*c*cos(lambda)^2)/30 + (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, (CLa*c*sin(2*lambda)*L^2)/120 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1)*L)/30, 0, 0, 0,     (CLa*L*c*sin(2*lambda))/20 + (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12 - (CLa*L^2*c*cos(lambda)^2)/20,                                -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0]
[ 0, 0,                                                                           0,                                                                             0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                                0, 0]];
    % Aerodynamic damping matrix
    Ca = [[ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                               (pi*c^2*cos(lambda)*sin(lambda))/2 - (26*pi*L*c*cos(lambda))/35 - (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),      (7*pi*L*c^2*cos(lambda)^2)/20 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (11*pi*L^2*c*cos(lambda))/105 + (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0, 0, 0,                                  (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L) - (9*pi*L*c*cos(lambda))/35 - (pi*c^2*cos(lambda)*sin(lambda))/2,    (3*pi*L*c^2*cos(lambda)^2)/20 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (13*pi*L^2*c*cos(lambda))/210 - (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0]
        [ 0, 0,                                                                     - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (7*pi*L*c*e*cos(lambda)^2)/10,                                               -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,                                                                       (pi*L^2*c*e*cos(lambda)^2)/10 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0, 0, 0,                                                                         (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (3*pi*L*c*e*cos(lambda)^2)/10,                                             -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,                                                                       (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192 - (pi*L^2*c*e*cos(lambda)^2)/15, 0]
        [ 0, 0, (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (11*pi*L^2*c*cos(lambda))/105 - (pi*L*c^2*cos(lambda)*sin(lambda))/10, - (pi*L^2*c^2*cos(lambda)^2)/20 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192,                                     - (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/60 - (2*pi*L^2*c^2*e*sin(lambda)^2)/15))/L - (2*pi*L^3*c*cos(lambda))/105, 0, 0, 0,   (13*pi*L^2*c*cos(lambda))/210 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (pi*L*c^2*cos(lambda)*sin(lambda))/10, (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192 - (pi*L^2*c^2*cos(lambda)^2)/30, (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/240 - (pi*L^2*c^2*e*sin(lambda)^2)/30))/L + (pi*L^3*c*cos(lambda))/70 + (pi*L^2*c^2*cos(lambda)*sin(lambda))/60, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                (pi*c^2*cos(lambda)*sin(lambda))/2 - (9*pi*L*c*cos(lambda))/35 + (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),      (3*pi*L*c^2*cos(lambda)^2)/20 + (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (13*pi*L^2*c*cos(lambda))/210 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0, 0, 0,                               - (pi*c^2*cos(lambda)*sin(lambda))/2 - (26*pi*L*c*cos(lambda))/35 - (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),    (7*pi*L*c^2*cos(lambda)^2)/20 + (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (pi*L*c^2*cos(lambda)*sin(lambda))/10 - (11*pi*L^2*c*cos(lambda))/105 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L, 0]
        [ 0, 0,                                                                     - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (3*pi*L*c*e*cos(lambda)^2)/10,                                               -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,                                                                       (pi*L^2*c*e*cos(lambda)^2)/15 + (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0, 0, 0,                                                                         (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (7*pi*L*c*e*cos(lambda)^2)/10,                                             -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,                                                                     - (pi*L^2*c*e*cos(lambda)^2)/10 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0]
        [ 0, 0, (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (13*pi*L^2*c*cos(lambda))/210 + (pi*L*c^2*cos(lambda)*sin(lambda))/10,   (pi*L^2*c^2*cos(lambda)^2)/30 + (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/240 - (pi*L^2*c^2*e*sin(lambda)^2)/30))/L + (pi*L^3*c*cos(lambda))/70 - (pi*L^2*c^2*cos(lambda)*sin(lambda))/60, 0, 0, 0, - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (11*pi*L^2*c*cos(lambda))/105 - (pi*L*c^2*cos(lambda)*sin(lambda))/10, (pi*L^2*c^2*cos(lambda)^2)/20 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192,                                     - (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/60 - (2*pi*L^2*c^2*e*sin(lambda)^2)/15))/L - (2*pi*L^3*c*cos(lambda))/105, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]];
    % Give a name to this
    fb = [                                                               0
        0
        (sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2 + (CLb*L*c*cos(lambda))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        -(CLb*L^2*c*cos(lambda))/12
        0
        0
        0
        (CLb*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        (CLb*L^2*c*cos(lambda))/12
        0];
    % Give a name to this
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0];
    % Give a name to this
    Lb =CLb*L*c*cos(lambda);
    % Vector fa (specify better the name )
    fa = [                                                          0;
        0;
        (sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2 + (CL0*L*c*cos(lambda))/2;
        (L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2;
        -(CL0*L^2*c*cos(lambda))/12;
        0;
        0;
        0;
        (CL0*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2;
        (L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2;
        (CL0*L^2*c*cos(lambda))/12;
        0];
    % Give a name to this
    lp = -(2*CLa*L*c*cos(lambda)^3*(x1^2 + x1*x2 + x2^2))/3;
    %
    lb = 2*CLb*L*c*cos(lambda)*(x1/2 + x2/2);
    %
    lq = [ 0, 0, -2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(2*x1 + x2))/3, -(CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0, 0, 0, 2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(x1 + 2*x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0];
    %
    Sq = [                                 0
        0
        (L*m*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(L*d*m*cos(lambda)^2*(2*x1 + x2))/6
        -(L^2*m*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        (L*m*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(L*d*m*cos(lambda)^2*(x1 + 2*x2))/6
        (L^2*m*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    %
    fp = [                                     0
        0
        -(CLa*L*c*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(2*x1 + x2))/6
        (CLa*L^2*c*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        -(CLa*L*c*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(x1 + 2*x2))/6
        -(CLa*L^2*c*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    
    
    %% Left wing
elseif dx==-1
    % Aerodynamic stiffness matrix
    Ka = [[ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]
        [ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]
        [ 0, 0, (CLa*c*sin(2*lambda))/4 - (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),      (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2 - (7*CLa*L*c*cos(lambda)^2)/20,      (CLa*L*c*sin(2*lambda))/20 + (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0, 0, 0,   (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L) - (CLa*c*sin(2*lambda))/4,    (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2 - (3*CLa*L*c*cos(lambda)^2)/20,      (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10 - (CLa*L*c*sin(2*lambda))/20, 0]
        [ 0, 0,                                 (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                     (CLa*L*c*e*cos(lambda)^3)/3,                                 (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                                  -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                   (CLa*L*c*e*cos(lambda)^3)/6,                                -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,   (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10 - (CLa*L*c*sin(2*lambda))/20,   (CLa*L^2*c*cos(lambda)^2)/20 + (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12,                              -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0, 0, 0,     (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, (CLa*L^2*c*cos(lambda)^2)/30 - (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, (CLa*c*sin(2*lambda)*L^2)/120 + (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1)*L)/30, 0]
        [ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]
        [ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]
        [ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]
        [ 0, 0, (CLa*c*sin(2*lambda))/4 + (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),    - (3*CLa*L*c*cos(lambda)^2)/20 - (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,    - (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0, 0, 0, - (CLa*c*sin(2*lambda))/4 - (6*CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/(5*L),  - (7*CLa*L*c*cos(lambda)^2)/20 - (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,      (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, 0]
        [ 0, 0,                                 (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                     (CLa*L*c*e*cos(lambda)^3)/6,                                -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                                  -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                                                   (CLa*L*c*e*cos(lambda)^3)/3,                                 (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,   (CLa*L*c*sin(2*lambda))/20 + (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, - (CLa*L^2*c*cos(lambda)^2)/30 - (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, (CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/30 - (CLa*L^2*c*sin(2*lambda))/120, 0, 0, 0,   - (CLa*L*c*sin(2*lambda))/20 - (CLa*c*e*cos(lambda)*(cos(lambda)^2 - 1))/10, (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12 - (CLa*L^2*c*cos(lambda)^2)/20,                              -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0]
        [ 0, 0,                                                                           0,                                                                               0,                                                                              0, 0, 0, 0,                                                                             0,                                                                             0,                                                                              0, 0]];
    % Aerodynamic damping matrix
    Ca = [[ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                               (pi*c^2*cos(lambda)*sin(lambda))/2 - (26*pi*L*c*cos(lambda))/35 - (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),      (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (7*pi*L*c^2*cos(lambda)^2)/20,    (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (11*pi*L^2*c*cos(lambda))/105 + (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0, 0, 0,                                  (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L) - (9*pi*L*c*cos(lambda))/35 - (pi*c^2*cos(lambda)*sin(lambda))/2,    (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 - (3*pi*L*c^2*cos(lambda)^2)/20,    (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (13*pi*L^2*c*cos(lambda))/210 - (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0]
        [ 0, 0,                                                                       (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 + (7*pi*L*c*e*cos(lambda)^2)/10,                                               -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,                                                                       (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192 - (pi*L^2*c*e*cos(lambda)^2)/10, 0, 0, 0,                                                                         (3*pi*L*c*e*cos(lambda)^2)/10 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,                                             -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,                                                                       (pi*L^2*c*e*cos(lambda)^2)/15 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0]
        [ 0, 0, (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (11*pi*L^2*c*cos(lambda))/105 - (pi*L*c^2*cos(lambda)*sin(lambda))/10,   (pi*L^2*c^2*cos(lambda)^2)/20 + (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192,                                     - (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/60 - (2*pi*L^2*c^2*e*sin(lambda)^2)/15))/L - (2*pi*L^3*c*cos(lambda))/105, 0, 0, 0,   (13*pi*L^2*c*cos(lambda))/210 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L + (pi*L*c^2*cos(lambda)*sin(lambda))/10, (pi*L^2*c^2*cos(lambda)^2)/30 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/240 - (pi*L^2*c^2*e*sin(lambda)^2)/30))/L + (pi*L^3*c*cos(lambda))/70 + (pi*L^2*c^2*cos(lambda)*sin(lambda))/60, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]
        [ 0, 0,                                (pi*c^2*cos(lambda)*sin(lambda))/2 - (9*pi*L*c*cos(lambda))/35 + (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),    - (3*pi*L*c^2*cos(lambda)^2)/20 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (13*pi*L^2*c*cos(lambda))/210 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (pi*L*c^2*cos(lambda)*sin(lambda))/10, 0, 0, 0,                               - (pi*c^2*cos(lambda)*sin(lambda))/2 - (26*pi*L*c*cos(lambda))/35 - (3*c^2*pi*cos(lambda)*sin(lambda)^2*(c - 8*e))/(20*L),  - (7*pi*L*c^2*cos(lambda)^2)/20 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,    (pi*L*c^2*cos(lambda)*sin(lambda))/10 - (11*pi*L^2*c*cos(lambda))/105 - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L, 0]
        [ 0, 0,                                                                       (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32 + (3*pi*L*c*e*cos(lambda)^2)/10,                                               -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,                                                                     - (pi*L^2*c*e*cos(lambda)^2)/15 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0, 0, 0,                                                                         (7*pi*L*c*e*cos(lambda)^2)/10 - (c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/32,                                             -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,                                                                       (pi*L^2*c*e*cos(lambda)^2)/10 + (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, 0]
        [ 0, 0, (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (13*pi*L^2*c*cos(lambda))/210 + (pi*L*c^2*cos(lambda)*sin(lambda))/10, - (pi*L^2*c^2*cos(lambda)^2)/30 - (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192, (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/240 - (pi*L^2*c^2*e*sin(lambda)^2)/30))/L + (pi*L^3*c*cos(lambda))/70 - (pi*L^2*c^2*cos(lambda)*sin(lambda))/60, 0, 0, 0, - (cos(lambda)*((pi*L*c^3*sin(lambda)^2)/80 - (pi*L*c^2*e*sin(lambda)^2)/10))/L - (11*pi*L^2*c*cos(lambda))/105 - (pi*L*c^2*cos(lambda)*sin(lambda))/10, (L*c^2*pi*sin(2*lambda)*cos(lambda)*(c - 8*e))/192 - (pi*L^2*c^2*cos(lambda)^2)/20,                                     - (cos(lambda)*((pi*L^2*c^3*sin(lambda)^2)/60 - (2*pi*L^2*c^2*e*sin(lambda)^2)/15))/L - (2*pi*L^3*c*cos(lambda))/105, 0]
        [ 0, 0,                                                                                                                                                     0,                                                                                    0,                                                                                                                                                        0, 0, 0, 0,                                                                                                                                                       0,                                                                                  0,                                                                                                                                                        0, 0]];
    % Give a name to this
    fb = [                                                               0
        0
        (sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2 + (CLb*L*c*cos(lambda))/2
        -(L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        -(CLb*L^2*c*cos(lambda))/12
        0
        0
        0
        (CLb*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2
        -(L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        (CLb*L^2*c*cos(lambda))/12
        0];
    % Give a name to this
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, -(CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, -(CLa*L*c*cos(lambda)^2)/2, 0, 0];
    % Give a name to this
    Lb =CLb*L*c*cos(lambda);
    % Vector fa (specify better the name )
    fa = [                                                          0;
        0;
        (sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2 + (CL0*L*c*cos(lambda))/2;
        -(L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2;
        -(CL0*L^2*c*cos(lambda))/12;
        0;
        0;
        0;
        (CL0*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2;
        -(L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2;
        (CL0*L^2*c*cos(lambda))/12;
        0];
    %% The following matrices and vectors are WRONG for a left wing, they
    % will be updated later
    lp = -(2*CLa*L*c*cos(lambda)^3*(x1^2 + x1*x2 + x2^2))/3;
    %
    lb = 2*CLb*L*c*cos(lambda)*(x1/2 + x2/2);
    %%
    lq = [ 0, 0, -2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(2*x1 + x2))/3, -(CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0, 0, 0, 2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(x1 + 2*x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0];
    %
    Sq = [                                 0
        0
        (L*m*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(L*d*m*cos(lambda)^2*(2*x1 + x2))/6
        -(L^2*m*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        (L*m*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(L*d*m*cos(lambda)^2*(x1 + 2*x2))/6
        (L^2*m*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    %
    fp = [                                     0
        0
        -(CLa*L*c*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(2*x1 + x2))/6
        (CLa*L^2*c*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        -(CLa*L*c*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(x1 + 2*x2))/6
        -(CLa*L^2*c*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    
else
    warning('The element is neither right nor left. Static aero matrices not created')
end

%% Assign matrices to the element

el.Ka=sparse(Ka);       % Stiffness matrix
el.Ca=sparse(Ca);       % Damping matrix
el.fa=sparse(fa);       % (Give a name to this)
el.lp = sparse(lp);     % (Give a name to this)
el.lb = sparse(lb);     % (Give a name to this)
el.lq = sparse(lq);     % (Give a name to this)
el.Sq = sparse(Sq);     % (Give a name to this)
el.fp = sparse(fp);     % (Give a name to this)
el.fb = fb;             % (Give a name to this)
el.Lq = Lq;             % (Give a name to this)
el.Lb = Lb;             % (Give a name to this)

end
