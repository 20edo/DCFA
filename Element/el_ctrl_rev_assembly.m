% Function that gives the expression of the aerodynamic terms needed for
% the computation of the control reversal problem
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
function el = el_ctrl_rev_assembly(b,el,lambda,alpha,dx)

L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; % TROVARE UN MODO PER ESPRIMERE e

%% Assign aero coefficients
CLb = el.sc.CLb;
CMb = el.sc.CMb;
CLa = el.sc.CLa;

%% Compute the matrices needed for control reversal

if dx == 1
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
    
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, (CLa*L*c*cos(lambda)^2)/2, 0, 0];
    
    Lb =CLb*L*c*cos(lambda);
elseif dx == -1
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
    
    Lq = [ 0, 0, (CLa*c*sin(2*lambda))/2, -(CLa*L*c*cos(lambda)^2)/2, 0, 0, 0, 0, -(CLa*c*sin(2*lambda))/2, -(CLa*L*c*cos(lambda)^2)/2, 0, 0];
    
    Lb =CLb*L*c*cos(lambda);
    
end 

el.fb = sparse(fb);
el.Lq = sparse(Lq);

el.fb = fb;
el.Lq = Lq;
el.Lb = Lb;

end