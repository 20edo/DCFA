function el=el_build_ssh_profile(b,insx,indx)
% Builds the element between two nodes, given the beam property.
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
%% Find the coordinate in which the section must be calculated
x=(insx.x+indx.x)/2;

%% Define constant section variables
c=b.c(x);
t=b.t(x);
h=b.h(x);
A=c*t*3;
E=b.E(x,0,0);
G=b.G(x,0,0);
rho=b.rho(x,0,0);
%% Compute section properties
sc=sc_ssh_profile(c,h,A,t,E,G,rho);
%% Set element properties
el = el_init();
el.sc=sc;
el.L=abs(indx.x-insx.x);
el=el_mass_assembly(el);
el=el_stiff_assembly(el);
el.C=zeros(size(el.M));
% Medium forces = forces on the medium section
el.navier=sparse([el.sc.navier  zeros(size(el.sc.navier))
                  zeros(size(el.sc.navier))  zeros(size(el.sc.navier))])*el.load;
end
