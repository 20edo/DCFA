% This is a function to compute the assembled matrices of the model
% starting from the beam matrices 
% The incence matrix of each beam is calculated

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
function model = m_compute_matrices(model)

n_node = size(model.en,2); % number of nodes of the structure

% dof_node has inside how many beams arrive or depart from one node
dof_node = zeros(1,n_node);
for i=1:size(model.b,2)
    a = model.b(i).on; 
    b = model.b(i).en; 
    dof_node(1,a) = dof_node(1,a) + 1;
    dof_node(1,b) = dof_node(1,b) + 1;
end

% compute the number of degrees of freedom of the entire structure
dof  = sum([model.b(:).nel]+1)-sum((dof_node(1,:)-1)); 
dof = 6*dof; 

% inizilise the matrices
model.M = sparse(zeros(dof));
model.K = sparse(zeros(dof));
model.C = sparse(zeros(dof));
model.Ka = sparse(zeros(dof));
model.f=  sparse(zeros(dof,1));


% organisation of dof: all external nodes and then all internal ones
% [6x6] node1 
%       [6x6] node2 
%            [6x6] node3 
%                 internal beam1 
%                               internal beam2    
%                                              ....

% Take care of the node matrices
index=zeros(1,6);
for i=1:n_node
    index=6*(i-1)+1;
    index=(index:index+5);
    model.M(index,index)=model.M(index,index)+model.en(i).M;
    model.K(index,index)=model.K(index,index)+model.en(i).K;
    model.C(index,index)=model.C(index,index)+model.en(i).C;
    model.f(index)=model.f(index)+model.en(i).f;
end 


% k is how many dof has been already added to the big matrix
k = 6*size(dof_node,2);


for i=1:size(model.b,2)
    a = model.b(i).on; 
    b = model.b(i).en; 
    n = size(model.b(i).M,1);
    % compute the incidence matrix of each the beam 
    model.b(i).A = sparse(zeros(dof,n));
    model.b(i).A(6*(a-1)+1:6*a,1:6) = eye(6); 
    model.b(i).A(6*(b-1)+1:6*b,end-5:end) = eye(6);
    model.b(i).A(k+1:k+n-12,7:end-6) = eye(n-12); 
    k = k+n-12; % refresh the number of dof used
    beam=b_rotate(model.b(i));
    % compute the model matrices
    model.M = model.M + (model.b(i).A)*(beam.M)*(model.b(i).A)';
    model.K = model.K + (model.b(i).A)*(beam.K)*(model.b(i).A)';
    model.C = model.C + (model.b(i).A)*(beam.C)*(model.b(i).A)';
    model.Ka = model.Ka + (model.b(i).A)*(beam.Ka)*(model.b(i).A)';
    model.f = model.f + (model.b(i).A)*(beam.f);
end

%% Constraint 
% explicit the external constraint eliminating rows and columns 
n_dv = 0; % number of contraint already applied to the model
for i=1:length(model.en)
    for k=1:6
        if model.en(i).c(k) % remove row and column 6(i-1)+k
            index = 6*(i-1)+k-n_dv;
            model.M = model.M([1:index-1,index+1:end],[1:index-1,index+1:end]);
            model.K = model.K([1:index-1,index+1:end],[1:index-1,index+1:end]);
            model.C = model.C([1:index-1,index+1:end],[1:index-1,index+1:end]);
            model.Ka = model.Ka([1:index-1,index+1:end],[1:index-1,index+1:end]);
            model.f = model.f([1:index-1,index+1:end]);
            n_dv = n_dv + 1; 
        end
    end            
end