function b=build_matrices(b)
% Builds the matrices of a beam, given the elements and nodes.
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
b.M=zeros(6*(b.nel+1));
for i=1:b.nel
    b.M(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))=b.M(6*(i-1)+1:6*(i+1),6*(i-1)+1:6*(i+1))+b.el(i).M.M;
end