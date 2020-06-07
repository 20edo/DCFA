% Builds the unsteady aero matrices of a beam.
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
% b->           Beam
% lambda->      Sweep angle
% dx->          1 for "right wings, -1 for "left" wings
% w_frac_vinf-> reduced frequency * chord (i.e. frequency/v_infinity)
function b=b_dynamic_aero_assembly(b,lambda,dx,k)
for j=1:b.nel
    %% Calculate element matrices
    b.el(j) = el_dynamic_aero_assembly(b,b.el(j),lambda,dx,k);
    %% Update beam matrices
    b.Ham(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))=b.Ham(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))+b.el(j).Ham;
    b.Hamb(6*(j-1)+1:6*(j+1),1)=b.Hamb(6*(j-1)+1:6*(j+1),1)+b.el(j).Hamb;
    b.Hbam(1,6*(j-1)+1:6*(j+1))=b.Hbam(1,6*(j-1)+1:6*(j+1))+b.el(j).Hbam;
    b.Hbb = b.Hbb + b.el(j).Hbb;
end



