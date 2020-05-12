structure = L_shaped_structure; 

n_node = 0; 
for i=1:size(structure.b,2)
    a = structure.b(i).on; 
    b = structure.b(i).en; 
    n_node = max([a,b,n_node]); 
end
dof_node = zeros(1,n_node);
for i=1:size(structure.b,2)
    a = structure.b(i).on; 
    b = structure.b(i).en; 
    dof_node(1,a) = dof_node(1,a) + 1;
    dof_node(1,b) = dof_node(1,b) + 1;
end

dof  = sum([structure.b(:).nel]+1)-sum((dof_node(1,:)-1)); 
dof = 6*dof; 
structure.M = zeros(dof);
structure.K = zeros(dof);
structure.C = zeros(dof);

% [6x6] node1 
%       [6x6] node2 
%            [6x6] node3 
%                 beam1 
%                      beam2    
%
k = 6*size(dof_node,2);
A = []; 
for i=1:size(structure.b,2)
    a = structure.b(i).on; 
    b = structure.b(i).en; 
    n = size(structure.b(i).M,1);
    structure.b(i).A = zeros(dof,n);
    structure.b(i).A(6*(a-1)+1:6*a,1:6) = eye(6); 
    structure.b(i).A(6*(b-1)+1:6*b,end-5:end) = eye(6);
    structure.b(i).A(k+1:k+n-12,7:end-6) = eye(n-12); 
    temp_Att = structure.b(i).A;
    temp_M = structure.b(i).M;
    k = k+n-12;
    structure.M = structure.M + (structure.b(i).A)*(structure.b(i).M)*(structure.b(i).A)';
    structure.K = structure.K + (structure.b(i).A)*(structure.b(i).K)*(structure.b(i).A)';
    structure.C = structure.C + (structure.b(i).A)*(structure.b(i).C)*(structure.b(i).A)';
end

















