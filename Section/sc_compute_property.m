function sc=sc_compute_property(sc)

% Calculate the elastic and insertial propreties of a section
% The function handles must be defined to accept arrays and return arrays
% with the same dimension (trick: add y.*0+z.*0).
% Torsion properties are not computed because they are very expensive, to
% compute torsion properties see sc_torsion
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
sc=sc_inertia(sc);
sc=sc_elastic(sc);

end