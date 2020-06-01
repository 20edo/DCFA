function el = el_ctrl_rev_assembly(b,el,lambda,alpha)

L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; % TROVARE UN MODO PER ESPRIMERE e

%% Compute CL/b
angle_b = [deg2rad(2), deg2rad(4), deg2rad(6)];
CL = [0.2232, 0.4440, 0.6586];
temp = polyfit(angle_b, CL, 1);
CLb = temp(1);
clear temp

%% Compute CM/b
angle_b = [deg2rad(2), deg2rad(4), deg2rad(6)];
CM = [-0.0309, -0.0613, -0.0898];
temp = polyfit(angle_b, CM, 1);
CMb = temp(1);
clear temp

%% 
CLa = 2*pi; % CL derived alpha


%% Compute the matrices needed for control reversal 

fb = [ [0];
    [0];
    [(sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2 + (CLb*L*c*cos(lambda))/2];
    [(L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2];
    [-(CLb*L^2*c*cos(lambda))/12];
    [0];
    [0];
    [0];
    [(CLb*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMb*c^2 + CLb*e*c))/2];
    [(L*cos(lambda)^2*(CMb*c^2 + CLb*e*c))/2];
    [(CLb*L^2*c*cos(lambda))/12];
    [0] ];

Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0];

Lb = CLb*L*c; 


el.fb = sparse(fb);
el.Lq = sparse(Lq);

el.fb = fb;
el.Lq = Lq;
el.Lb = Lb; 

end