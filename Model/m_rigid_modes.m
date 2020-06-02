% This function returns the eigenvectors of the rigid body modes
% The model MUST be free-free (no external constraints)
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

 

function V = m_rigid_modes(m)
V  = zeros(size(m.K,1),6); 
% Trasl along x
V(1:6:end,1) = 1;  
V(:,1) = V(:,1)/norm(V(:,1)); 
% Trasl along y
V(2:6:end,2) = 1;  
V(:,2) = V(:,2)/norm(V(:,2));
% Trasl along z
V(3:6:end,3) = 1;  
V(:,3) = V(:,3)/norm(V(:,3));
% Rotations are defined wrt [0,0,0]'
alpha=2;
%% External nodes
% Rotation around x
for i=1:length(m.en)
   V(6*(i-1)+1:6*i,4)=[rotx(alpha)*m.en(i).x-m.en(i).x; deg2rad(alpha); 0; 0];
end

 

%  Rotation around y
for i=1:length(m.en)
   V(6*(i-1)+1:6*i,5)=[roty(alpha)*m.en(i).x-m.en(i).x; 0; deg2rad(alpha); 0];
end

 

%  Rotation around z
for i=1:length(m.en)
   V(6*(i-1)+1:6*i,6)=[rotz(alpha)*m.en(i).x-m.en(i).x; 0; 0; deg2rad(alpha)];
end

 


clear i
%% Internal nodes
k = 6*length(m.en);
for j=1:length(m.b)
    for i=2:length(m.b(j).in)-1
        a = i-1; 
        pos = m.b(j).o + m.b(j).vx*m.b(j).in(i).x;
        V(k+6*(a-1)+1:k+6*a,4)=[rotx(alpha)*pos-pos; deg2rad(alpha); 0; 0];
        V(k+6*(a-1)+1:k+6*a,5)=[roty(alpha)*pos-pos; 0; deg2rad(alpha); 0];
        V(k+6*(a-1)+1:k+6*a,6)=[rotz(alpha)*pos-pos; 0; 0; deg2rad(alpha)];
    end
        k = k + 6*(length(m.b(j).in)-2);
end
V(:,4) = V(:,4)/norm(V(:,4));        
V(:,5) = V(:,5)/norm(V(:,5));     
V(:,6) = V(:,6)/norm(V(:,6));

end 
    