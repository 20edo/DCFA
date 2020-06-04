function b=build_matrices(b)
% Builds the matrices of a beam, given the elements and nodes.
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
b.M=sparse(zeros(6*(b.nel+1)));
b.K=sparse(zeros(6*(b.nel+1)));
b.C=sparse(zeros(6*(b.nel+1)));
b.Ka=sparse(zeros(6*(b.nel+1)));
b.Ca=sparse(zeros(6*(b.nel+1)));
b.f=sparse(zeros(6*(b.nel+1),1));
b.fa=sparse(zeros(6*(b.nel+1),1));
b.fb=sparse(zeros(6*(b.nel+1),1));
b.Lq=sparse(zeros(1,6*(b.nel+1)));
b.Lb = 0; 
b.Jx = 0; 
b.lp = 0; 
b.lb = 0;
b.lq=sparse(zeros(1,6*(b.nel+1)));
b.Sq=sparse(zeros(6*(b.nel+1),1));
b.fp=sparse(zeros(6*(b.nel+1),1));


for i=1:b.nel
    b.M(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))=b.M(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))+b.el(i).M;
    b.K(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))=b.K(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))+b.el(i).K;
    b.C(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))=b.C(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))+b.el(i).C;   
end