function el = el_consistent_assembly(b,el,i,lambda)

L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; 
d = el.sc.ycg;
m = el.sc.m; 

pos1 = b.o + b.vx*b.in(i).x;
pos2 = b.o + b.vx*b.in(i+1).x;
x1 = pos1(2);
x2 = pos2(2);

CLa = 2*pi; 
 
%% Compute CM/b
angle_b = [deg2rad(2), deg2rad(4), deg2rad(6)];
CM = [-0.0309, -0.0613, -0.0898]; % 70% alpha zero
%CM = [-0.0275, -0.0546, -0.0802]; % 85% alpha zero
temp = polyfit(angle_b, CM, 1);
CMb = temp(1); 
CMb = -1.7; 
% CMb = -0.75;
clear temp
%% Compute CL/b
angle_b = [deg2rad(2), deg2rad(4), deg2rad(6)];
% CL = [0.2232, 0.4440, 0.6586]; % 70% alpha zero
CL = [0.1544, 0.3073, 0.4534]; % 85% alpha zero
temp = polyfit(angle_b, CL, 1);
CLb = temp(1);
CLb = 0.7;
% CLb = 0.3;
clear temp


%%
Jx = (2*L*m*cos(lambda)^3*(x1^2 + x1*x2 + x2^2))/3;
%%
lp = -(2*CLa*L*c*cos(lambda)^3*(x1^2 + x1*x2 + x2^2))/3;
%%
lb = 2*CLb*L*c*cos(lambda)*(x1/2 + x2/2);
%%
lq = [ 0, 0, -2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(2*x1 + x2))/3, -(CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0, 0, 0, 2*CLa*c*sin(lambda)*(x1/2 + x2/2)*(sin(lambda)^2 - 1), (CLa*L*c*cos(lambda)^3*(x1 + 2*x2))/3, (CLa*L*c*sin(lambda)*(x1/2 - x2/2)*(sin(lambda)^2 - 1))/3, 0];
%%
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
%%
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

el.Jx = sparse(Jx); 
el.lp = sparse(lp);
el.lb = sparse(lb);
el.lq = sparse(lq);
el.Sq = sparse(Sq);
el.fp = sparse(fp);

