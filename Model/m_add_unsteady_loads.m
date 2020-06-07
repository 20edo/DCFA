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
function m = m_add_unsteady_loads(m,v_inf,k)

for i = 1:length(m.b)
    if m.b(i).ssh % Check if the beam is a profile
        vz = cross(m.b(i).vx,m.b(i).vy);
        % Calculate sweep angle
        lambda = -(acos(dot(m.b(i).vx,v_inf/norm(v_inf)))-pi/2);
        % Check if it is a dx profile or sx profile
        dx = sign(m.b(i).vx(2)); 
        % in the case of the ruddder, the left or right behaviour depends
        % on the direction of the wind 
        if dx == 0
            dx = - sign(v_inf(2));              
        end
        if dx == 0
            disp(m.b(i).name);
            dx = 1;
        end
        m.b(i)=b_dynamic_aero_assembly(m.b(i),lambda,dx,k);
    end
end


n_node = size(m.en,2); % number of nodes of the structure

% dof_node has inside how many beams arrive or depart from one node
dof_node = zeros(1,n_node);
for i=1:size(m.b,2)
    a = m.b(i).on; 
    b = m.b(i).en; 
    dof_node(1,a) = dof_node(1,a) + 1;
    dof_node(1,b) = dof_node(1,b) + 1;
end

% compute the number of degrees of freedom of the entire structure
dof  = sum([m.b(:).nel]+1)-sum((dof_node(1,:)-1)); 
dof = 6*dof; 

% inizilise the matrices
m.Ham = sparse(zeros(dof));
m.Hamb=  sparse(zeros(dof,1));
m.Hbam = sparse(zeros(1,dof));
m.Hbb = 0; 


% organisation of dof: all external nodes and then all internal ones
% [6x6] node1 
%       [6x6] node2 
%            [6x6] node3 
%                 internal beam1 
%                               internal beam2    
%                                              ....

for i=1:size(m.b,2)
    beam=b_rotate(m.b(i));
    % compute the model matrices
    m.Ham = m.Ham + (m.b(i).A)*(beam.Ham)*(m.b(i).A)';
    m.Hamb = m.Hamb + (m.b(i).A)*(beam.Hamb);
    m.Hbam = m.Hbam + (beam.Hbam)*(m.b(i).A)';
    m.Hbb = m.Hbb + beam.Hbb; 
end

%% Constraint 
% explicit the external constraint eliminating rows and columns 
n_dv = 0; % number of contraint already applied to the model
for i=1:length(m.en)
    for k=1:6
        if m.en(i).c(k) % remove row and column 6(i-1)+k
            index = 6*(i-1)+k-n_dv;
            m.Ham = m.Ham([1:index-1,index+1:end],[1:index-1,index+1:end]);
            m.Hamb = m.Hamb([1:index-1,index+1:end]);
            m.Hbam = m.Hbam([1:index-1,index+1:end]);
            n_dv = n_dv + 1; 
        end
    end            
end