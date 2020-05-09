function sc=sc_constant_p_rect(a,b,E,G,rho)

% Returns a rectangular section of sides a(along y) and b (along z) and constant properties
% Improve this function: write the exact values instad of integrating
% numerically
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
%% Data structures definition:

sc=sc_init();
sc.Zmax=@(y) b/2+y.*0;
sc.Zmin=@(y) -b/2+y.*0;
sc.Ymax=a/2;
sc.Ymin=-a/2;
sc.pol=false;               % True if the section is defined in polar coordinates
sc.E=@(y,z) E+0.*y+0.*z;
sc.G=@(y,z) G+0.*y+0.*z;
sc.rho=@(y,z) rho+0.*y+0.*z;
sc=sc_compute_property(sc);

end
