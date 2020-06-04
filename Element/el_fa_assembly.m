% Function that creates the vector fa of the static aero loads of the
% structure. It is divided in two parts: 
% dx = 1  -> right wing 
% dx = -1 -> left wing 
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
function el = el_fa_assembly(b,el,lambda,alpha,dx)

L=el.L;
c = el.sc.Ymax - el.sc.Ymin ;  % CORDA
e = -el.sc.yo - 1/4*c; % TROVARE UN MODO PER ESPRIMERE e

CLa=el.sc.CLa;
CL0 = alpha*CLa; % CL at zero incidence
CMac = el.sc.CMac; % Moment coefficient wrt aerodynamic center

if dx == 1
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
elseif dx == -1
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
end

el.fa=sparse(fa);
el.fa=fa;

end
                                                            
                                                            
                                                            
                                                            
                                                          
                                                            
                                                            