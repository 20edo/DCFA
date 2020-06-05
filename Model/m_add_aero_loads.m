% Compute the static aerodynamic matrices
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
    if m.b(i).ssh % Check if the beam is a profile
        vz = cross(m.b(i).vx,m.b(i).vy);
        % Calculate alpha
        alpha = sign(m.b(i).vx(2))*(atan2(m.b(i).vy'*v_inf/norm(v_inf),vz'*v_inf/norm(v_inf))-pi/2);
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
        m.b(i)=b_static_aero_assembly(m.b(i),lambda,alpha,dx);
    end
end
m = m_compute_matrices(m);



