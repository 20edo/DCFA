function b=b_static_aero_assembly(b,lambda,alpha,dx)
% Builds the static aero matrices of a beam.
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

for j=1:b.nel
    %% Calculate element matrices
    b.el(j) = el_static_aero_assembly(b,b.el(j),lambda,j,alpha,dx);
    %% Update beam matrices
    b.Ka(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))=b.Ka(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))+b.el(j).Ka;
    b.Ca(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))=b.Ca(6*(j-1)+1:6*(j+1),6*(j-1)+1:6*(j+1))+b.el(j).Ca;
    b.fa(6*(j-1)+1:6*(j+1),1)=b.fa(6*(j-1)+1:6*(j+1),1)+b.el(j).fa;
    b.fb(6*(j-1)+1:6*(j+1),1)=b.fb(6*(j-1)+1:6*(j+1),1)+b.el(j).fb;
    b.Lq(1,6*(j-1)+1:6*(j+1))=b.Lq(1,6*(j-1)+1:6*(j+1))+b.el(j).Lq;
    b.Lb = b.Lb + b.el(j).Lb;
    b.lp = b.lp + b.el(j).lp;
    b.lb = b.lb + b.el(j).lb;
    b.lq(1,6*(j-1)+1:6*(j+1))=b.lq(1,6*(j-1)+1:6*(j+1))+b.el(j).lq;
    b.Sq(6*(j-1)+1:6*(j+1),1)=b.Sq(6*(j-1)+1:6*(j+1),1)+b.el(j).Sq;
    b.fp(6*(j-1)+1:6*(j+1),1)=b.fp(6*(j-1)+1:6*(j+1),1)+b.el(j).fp;
end



