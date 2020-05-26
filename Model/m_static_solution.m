function m=m_static_solution(m)
% Calculate the static solution of the problem

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%

u=m.K\m.f;

%% Put the displacements in the appropriate field

for j=1:length(m.en)
    for k=1:6
        if m.en(j).c(k) % remove row and column 6(i-1)+k
            index = 6*(j-1)+k;
            u = [u(1:index-1,:);zeros(1,size(u,2));u(index:end,:)]; 
        end
    end            
end

for i = 1:length(m.b)
    u_beam = transpose(m.b(i).A)*u; 
    for k = 1:m.b(i).nel+1
%        model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),end);
%        if we would like to reconstruct the solution over time
       m.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end

end

