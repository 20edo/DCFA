% Compute the aerodynamic stiffness matrix of the model (so far)
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
function m=m_add_aero_loads(m,v_inf)

for i = 1:length(m.b)
    if m.b(i).ssh
        vz = cross(m.b(i).vx,m.b(i).vy);
        % Calculate alpha
        alpha = sign(m.b(i).vx(2))*(atan2(m.b(i).vy'*v_inf/norm(v_inf),vz'*v_inf/norm(v_inf))-pi/2);
        % Calculate sweep angle
        lambda = -(acos(dot(m.b(i).vx,v_inf/norm(v_inf)))-pi/2);
        % Check if the beam is a profile
        for j=1:m.b(i).nel
            m.b(i).el(j) = el_Ka_assembly(m.b(i),m.b(i).el(j),lambda,alpha);
            % Calculate K_aero
        end
        for k=1:m.b(i).nel % Create Ka matrix for the beam
            m.b(i).Ka(6*(k-1)+1:6*(k+1),6*(k-1)+1:6*(k+1))=m.b(i).Ka(6*(k-1)+1:6*(k+1),6*(k-1)+1:6*(k+1))+m.b(i).el(k).Ka;
        end
    end
end
m = m_compute_matrices(m);



