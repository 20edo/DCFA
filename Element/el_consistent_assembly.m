function el = el_consistent_assembly(b,el,i,lambda)

%% Define geometrical variables
L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; 
d = el.sc.ycg;
m = el.sc.m; 

pos1 = b.o + b.vx*b.in(i).x;
pos2 = b.o + b.vx*b.in(i+1).x;
x1 = pos1(2);
x2 = pos2(2);

 
%% Recover aero coefficients

CLb=el.sc.CLb;
CMb=el.sc.CMb;
CLa = el.sc.CLa;

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

