function el = el_Ka_assembly(b,el,lambda,alpha)

L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; % TROVARE UN MODO PER ESPRIMERE e


CL0 = alpha*2*pi; % CL at zero incidence
CMac = -pi/2*(-alpha); % CM wrt the aerodynamic center

fa = [[0];
    [0];
    [(sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2 + (CL0*L*c*cos(lambda))/2];
    [(L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2];
    [-(CL0*L^2*c*cos(lambda))/12];
    [0];
    [0];
    [0];
    [(CL0*L*c*cos(lambda))/2 - (sin(2*lambda)*(CMac*c^2 + CL0*e*c))/2];
    [(L*cos(lambda)^2*(CMac*c^2 + CL0*e*c))/2];
    [(CL0*L^2*c*cos(lambda))/12];
    [0]];

el.fa=sparse(fa);
el.fa=fa;

end
                                                            
                                                            
                                                            
                                                            
                                                          
                                                            
                                                            