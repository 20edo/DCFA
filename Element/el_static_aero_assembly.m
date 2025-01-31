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
    Ka = [[ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0, (CLa*c*cos(lambda)*sin(lambda)*(5*L + 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(7*L + 10*e*sin(lambda)))/20,       (CLa*c*sin(2*lambda)*(L - e*sin(lambda)))/20, 0, 0, 0, -(CLa*c*cos(lambda)*sin(lambda)*(5*L + 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(3*L + 10*e*sin(lambda)))/20,     -(CLa*c*sin(2*lambda)*(L + e*sin(lambda)))/20, 0]
        [ 0, 0,                     (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/3,     (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                     -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/6,   -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,         -(CLa*c*cos(lambda)*sin(lambda)*(L + e*sin(lambda)))/10, -(CLa*L*c*(sin(lambda)^2 - 1)*(3*L - 5*e*sin(lambda)))/60,  -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0, 0, 0,           (CLa*c*cos(lambda)*sin(lambda)*(L + e*sin(lambda)))/10, -(CLa*L*c*(sin(lambda)^2 - 1)*(2*L + 5*e*sin(lambda)))/60, (CLa*L*c*sin(2*lambda)*(L - 2*e*sin(lambda)))/120, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0, (CLa*c*cos(lambda)*sin(lambda)*(5*L - 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(3*L - 10*e*sin(lambda)))/20,      -(CLa*c*sin(2*lambda)*(L - e*sin(lambda)))/20, 0, 0, 0, -(CLa*c*cos(lambda)*sin(lambda)*(5*L - 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(7*L - 10*e*sin(lambda)))/20,      (CLa*c*sin(2*lambda)*(L + e*sin(lambda)))/20, 0]
        [ 0, 0,                     (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/6,    -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                     -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/3,    (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,          (CLa*c*cos(lambda)*sin(lambda)*(L - e*sin(lambda)))/10,  (CLa*L*c*(sin(lambda)^2 - 1)*(2*L - 5*e*sin(lambda)))/60, -(CLa*L*c*sin(2*lambda)*(L + 2*e*sin(lambda)))/120, 0, 0, 0,          -(CLa*c*cos(lambda)*sin(lambda)*(L - e*sin(lambda)))/10,  (CLa*L*c*(sin(lambda)^2 - 1)*(3*L + 5*e*sin(lambda)))/60, -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]];
    % Aerodynamic damping matrix
    Ca = [[ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0, -(L*cos(lambda)*((52*pi*c)/35 - (c*pi*(35*c^2*sin(lambda)^2 + 280*c*sin(lambda) - 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 56*L*cos(lambda)))/160, (L*c*pi*cos(lambda)*(176*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, 0, 0, 0, -(L*cos(lambda)*((18*pi*c)/35 + (c*pi*(280*c*sin(lambda) - 35*c^2*sin(lambda)^2 + 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 24*L*cos(lambda)))/160, -(L*c*pi*cos(lambda)*(104*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, 0]
        [ 0, 0,             (L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L + (7*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,        -(L*c*pi*cos(lambda)*(- 5*sin(2*lambda)*c^2 + 40*e*sin(2*lambda)*c + 96*L*e*cos(lambda)))/960, 0, 0, 0,            -(L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L - (3*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,          (L*c*pi*cos(lambda)*(- 5*sin(2*lambda)*c^2 + 40*e*sin(2*lambda)*c + 64*L*e*cos(lambda)))/960, 0]
        [ 0, 0,                  (L*c*pi*cos(lambda)*(176*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680,  (L*c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 48*L*cos(lambda)))/960,                                                                        -(2*pi*L^3*c*cos(lambda))/105, 0, 0, 0,                  (L*c*pi*cos(lambda)*(104*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680,  (L*c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 32*L*cos(lambda)))/960,    (L^2*c*pi*cos(lambda)*(48*L - 7*c^2*sin(lambda)^2 + 56*c*sin(lambda) + 56*c*e*sin(lambda)^2))/3360, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0, -(L*cos(lambda)*((18*pi*c)/35 - (c*pi*(280*c*sin(lambda) - 35*c^2*sin(lambda)^2 + 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 24*L*cos(lambda)))/160, (L*c*pi*cos(lambda)*(104*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, 0, 0, 0, -(L*cos(lambda)*((52*pi*c)/35 + (c*pi*(35*c^2*sin(lambda)^2 + 280*c*sin(lambda) - 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 56*L*cos(lambda)))/160, -(L*c*pi*cos(lambda)*(176*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, 0]
        [ 0, 0,             (L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L + (3*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,          -(L*c*pi*cos(lambda)*(5*sin(2*lambda)*c^2 - 40*e*sin(2*lambda)*c + 64*L*e*cos(lambda)))/960, 0, 0, 0,            -(L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L - (7*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,            (L*c*pi*cos(lambda)*(5*sin(2*lambda)*c^2 - 40*e*sin(2*lambda)*c + 96*L*e*cos(lambda)))/960, 0]
        [ 0, 0,                 -(L*c*pi*cos(lambda)*(104*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, -(L*c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 32*L*cos(lambda)))/960,   (L^2*c*pi*cos(lambda)*(48*L + 7*c^2*sin(lambda)^2 - 56*c*sin(lambda) - 56*c*e*sin(lambda)^2))/3360, 0, 0, 0,                 -(L*c*pi*cos(lambda)*(176*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, -(L*c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 48*L*cos(lambda)))/960,                                                                         -(2*pi*L^3*c*cos(lambda))/105, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]];
    % Give a name to this
    fb = [                                                        0
        0
        -(L*((sin(2*lambda)*(CMb*c^2 + CLb*e*c))/L + CLb*c*cos(lambda)))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        (CLb*L^2*c*cos(lambda))/12
        0
        0
        0
        (L*((sin(2*lambda)*(CMb*c^2 + CLb*e*c))/L - CLb*c*cos(lambda)))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        -(CLb*L^2*c*cos(lambda))/12
        0];
    % Give a name to this
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0];
    
    % Give a name to this
    Lb = CLb*L*c*cos(lambda);
    
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
    lq = -[ 0, 0, 2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), -(CLa*L*c*cos(lambda)^3*(2*x1 + x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0, 0, 0, -2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(x1 + 2*x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0];
    Sq = [                             0
        0
        -(L*m*cos(lambda)^2*(7*x1 + 3*x2))/20
        (L*d*m*cos(lambda)^2*(2*x1 + x2))/6
        (L^2*m*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        -(L*m*cos(lambda)^2*(3*x1 + 7*x2))/20
        (L*d*m*cos(lambda)^2*(x1 + 2*x2))/6
        -(L^2*m*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    %
    fp = [                                 0
        0
        (CLa*L*c*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(2*x1 + x2))/6
        -(CLa*L^2*c*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        (CLa*L*c*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(x1 + 2*x2))/6
        (CLa*L^2*c*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    
    
    %% Left wing
elseif dx==-1
    % Aerodynamic stiffness matrix
    Ka = [[ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0, (CLa*c*cos(lambda)*sin(lambda)*(5*L + 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(7*L + 10*e*sin(lambda)))/20,       (CLa*c*sin(2*lambda)*(L - e*sin(lambda)))/20, 0, 0, 0, -(CLa*c*cos(lambda)*sin(lambda)*(5*L + 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(3*L + 10*e*sin(lambda)))/20,     -(CLa*c*sin(2*lambda)*(L + e*sin(lambda)))/20, 0]
        [ 0, 0,                     (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/3,     (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                     -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/6,   -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,         -(CLa*c*cos(lambda)*sin(lambda)*(L + e*sin(lambda)))/10, -(CLa*L*c*(sin(lambda)^2 - 1)*(3*L - 5*e*sin(lambda)))/60,  -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0, 0, 0,           (CLa*c*cos(lambda)*sin(lambda)*(L + e*sin(lambda)))/10, -(CLa*L*c*(sin(lambda)^2 - 1)*(2*L + 5*e*sin(lambda)))/60, (CLa*L*c*sin(2*lambda)*(L - 2*e*sin(lambda)))/120, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]
        [ 0, 0, (CLa*c*cos(lambda)*sin(lambda)*(5*L - 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(3*L - 10*e*sin(lambda)))/20,      -(CLa*c*sin(2*lambda)*(L - e*sin(lambda)))/20, 0, 0, 0, -(CLa*c*cos(lambda)*sin(lambda)*(5*L - 12*e*sin(lambda)))/(10*L),   (CLa*c*(sin(lambda)^2 - 1)*(7*L - 10*e*sin(lambda)))/20,      (CLa*c*sin(2*lambda)*(L + e*sin(lambda)))/20, 0]
        [ 0, 0,                     (CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/6,    -(CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0, 0, 0,                     -(CLa*c*e*sin(lambda)*(sin(lambda)^2 - 1))/2,                               (CLa*L*c*e*cos(lambda)^3)/3,    (CLa*L*c*e*sin(lambda)*(sin(lambda)^2 - 1))/12, 0]
        [ 0, 0,          (CLa*c*cos(lambda)*sin(lambda)*(L - e*sin(lambda)))/10,  (CLa*L*c*(sin(lambda)^2 - 1)*(2*L - 5*e*sin(lambda)))/60, -(CLa*L*c*sin(2*lambda)*(L + 2*e*sin(lambda)))/120, 0, 0, 0,          -(CLa*c*cos(lambda)*sin(lambda)*(L - e*sin(lambda)))/10,  (CLa*L*c*(sin(lambda)^2 - 1)*(3*L + 5*e*sin(lambda)))/60, -(2*CLa*L*c*e*cos(lambda)*(cos(lambda)^2 - 1))/15, 0]
        [ 0, 0,                                                               0,                                                         0,                                                  0, 0, 0, 0,                                                                0,                                                         0,                                                 0, 0]];
    
    Ka([4,10],:) = -Ka([4,10],:);
    Ka(:,[4,10]) = -Ka(:,[4,10]);
    
    % Aerodynamic damping matrix
    Ca = [[ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0, -(L*cos(lambda)*((52*pi*c)/35 - (c*pi*(35*c^2*sin(lambda)^2 + 280*c*sin(lambda) - 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 56*L*cos(lambda)))/160, (L*c*pi*cos(lambda)*(176*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, 0, 0, 0, -(L*cos(lambda)*((18*pi*c)/35 + (c*pi*(280*c*sin(lambda) - 35*c^2*sin(lambda)^2 + 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 24*L*cos(lambda)))/160, -(L*c*pi*cos(lambda)*(104*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, 0]
        [ 0, 0,             (L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L + (7*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,        -(L*c*pi*cos(lambda)*(- 5*sin(2*lambda)*c^2 + 40*e*sin(2*lambda)*c + 96*L*e*cos(lambda)))/960, 0, 0, 0,            -(L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L - (3*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,          (L*c*pi*cos(lambda)*(- 5*sin(2*lambda)*c^2 + 40*e*sin(2*lambda)*c + 64*L*e*cos(lambda)))/960, 0]
        [ 0, 0,                  (L*c*pi*cos(lambda)*(176*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680,  (L*c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 48*L*cos(lambda)))/960,                                                                        -(2*pi*L^3*c*cos(lambda))/105, 0, 0, 0,                  (L*c*pi*cos(lambda)*(104*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680,  (L*c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 32*L*cos(lambda)))/960,    (L^2*c*pi*cos(lambda)*(48*L - 7*c^2*sin(lambda)^2 + 56*c*sin(lambda) + 56*c*e*sin(lambda)^2))/3360, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]
        [ 0, 0, -(L*cos(lambda)*((18*pi*c)/35 - (c*pi*(280*c*sin(lambda) - 35*c^2*sin(lambda)^2 + 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 24*L*cos(lambda)))/160, (L*c*pi*cos(lambda)*(104*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, 0, 0, 0, -(L*cos(lambda)*((52*pi*c)/35 + (c*pi*(35*c^2*sin(lambda)^2 + 280*c*sin(lambda) - 280*c*e*sin(lambda)^2))/(280*L)))/2,   -(c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 56*L*cos(lambda)))/160, -(L*c*pi*cos(lambda)*(176*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, 0]
        [ 0, 0,             (L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L + (3*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/48,          -(L*c*pi*cos(lambda)*(5*sin(2*lambda)*c^2 - 40*e*sin(2*lambda)*c + 64*L*e*cos(lambda)))/960, 0, 0, 0,            -(L*cos(lambda)*(((pi*sin(2*lambda)*c^3)/16 - (e*pi*sin(2*lambda)*c^2)/2)/L - (7*pi*c*e*cos(lambda))/5))/2,                                                  -(L*c^2*pi*cos(lambda)^3*(c - 8*e))/24,            (L*c*pi*cos(lambda)*(5*sin(2*lambda)*c^2 - 40*e*sin(2*lambda)*c + 96*L*e*cos(lambda)))/960, 0]
        [ 0, 0,                 -(L*c*pi*cos(lambda)*(104*L + 21*c^2*sin(lambda)^2 - 168*c*sin(lambda) - 168*c*e*sin(lambda)^2))/1680, -(L*c^2*pi*cos(lambda)*(5*c*sin(2*lambda) - 40*e*sin(2*lambda) + 32*L*cos(lambda)))/960,   (L^2*c*pi*cos(lambda)*(48*L + 7*c^2*sin(lambda)^2 - 56*c*sin(lambda) - 56*c*e*sin(lambda)^2))/3360, 0, 0, 0,                 -(L*c*pi*cos(lambda)*(176*L - 21*c^2*sin(lambda)^2 + 168*c*sin(lambda) + 168*c*e*sin(lambda)^2))/1680, -(L*c^2*pi*cos(lambda)*(40*e*sin(2*lambda) - 5*c*sin(2*lambda) + 48*L*cos(lambda)))/960,                                                                         -(2*pi*L^3*c*cos(lambda))/105, 0]
        [ 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                    0, 0, 0, 0,                                                                                                                     0,                                                                                       0,                                                                                                     0, 0]];
    
    Ca([4,10],:) = -Ca([4,10],:);
    Ca(:,[4,10]) = -Ca(:,[4,10]);
    % Give a name to this
    fb = [                                                        0
        0
        -(L*((sin(2*lambda)*(CMb*c^2 + CLb*e*c))/L + CLb*c*cos(lambda)))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        (CLb*L^2*c*cos(lambda))/12
        0
        0
        0
        (L*((sin(2*lambda)*(CMb*c^2 + CLb*e*c))/L - CLb*c*cos(lambda)))/2
        (L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2
        -(CLb*L^2*c*cos(lambda))/12
        0];
    
    
    fb([4,10]) = -fb([4,10]);
    % Give a name to this
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0];
    
    Lq([4,10]) = -Lq([4,10]);
    % Give a name to this
    Lb = CLb*L*c*cos(lambda);
    
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
    
    fa([4,10]) = -fa([4,10]);
    % Give a name to this
    lp = -(2*CLa*L*c*cos(lambda)^3*(x1^2 + x1*x2 + x2^2))/3;
    %
    lb = 2*CLb*L*c*cos(lambda)*(x1/2 + x2/2);
    %
    lq = -[ 0, 0, 2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), -(CLa*L*c*cos(lambda)^3*(2*x1 + x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0, 0, 0, -2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(x1 + 2*x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0];
    lq([4,10]) = -lq([4,10]);
    Sq = [                             0
        0
        -(L*m*cos(lambda)^2*(7*x1 + 3*x2))/20
        (L*d*m*cos(lambda)^2*(2*x1 + x2))/6
        (L^2*m*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        -(L*m*cos(lambda)^2*(3*x1 + 7*x2))/20
        (L*d*m*cos(lambda)^2*(x1 + 2*x2))/6
        -(L^2*m*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    Sq([4,10]) = -Sq([4,10]);
    %
    fp = [                                 0
        0
        (CLa*L*c*cos(lambda)^2*(7*x1 + 3*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(2*x1 + x2))/6
        -(CLa*L^2*c*cos(lambda)^2*(3*x1 + 2*x2))/60
        0
        0
        0
        (CLa*L*c*cos(lambda)^2*(3*x1 + 7*x2))/20
        -(CLa*L*c*e*cos(lambda)^2*(x1 + 2*x2))/6
        (CLa*L^2*c*cos(lambda)^2*(2*x1 + 3*x2))/60
        0];
    fp([4,10]) = -fp([4,10]);
    
else
    warning('The element is neither right nor left. Static aero matrices not created')
end

%% Assign matrices to the element
if isequal(b.name,'rudder')
    R = blkdiag(rotx(90),rotx(90),rotx(90),rotx(90));
    Ka = R*Ka*R'; 
end


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
