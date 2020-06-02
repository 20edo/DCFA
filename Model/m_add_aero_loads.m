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
  
        
        for j=1:m.b(i).nel
            % Calculate K_aero
            m.b(i).el(j) = el_Ka_assembly(m.b(i),m.b(i).el(j),lambda,alpha,dx);
            % Calculate f_aero
            m.b(i).el(j) = el_fa_assembly(m.b(i),m.b(i).el(j),lambda,alpha,dx);
            % Calculate matrices for control reversal 
            m.b(i).el(j) = el_ctrl_rev_assembly(m.b(i),m.b(i).el(j),lambda,alpha,dx);
            % Calculate matrices for consistend roll problems 
            m.b(i).el(j) = el_consistent_assembly(m.b(i),m.b(i).el(j),j,lambda);
        end
        for k=1:m.b(i).nel 
            % Create Ka matrix for the beam
            m.b(i).Ka(6*(k-1)+1:6*(k+1),6*(k-1)+1:6*(k+1))=m.b(i).Ka(6*(k-1)+1:6*(k+1),6*(k-1)+1:6*(k+1))+m.b(i).el(k).Ka;
            % Create the fa vector for the beam 
            m.b(i).fa(6*(k-1)+1:6*(k+1),1)=m.b(i).fa(6*(k-1)+1:6*(k+1),1)+m.b(i).el(k).fa;
            % Create the fb vector for the beam 
            m.b(i).fb(6*(k-1)+1:6*(k+1),1)=m.b(i).fb(6*(k-1)+1:6*(k+1),1)+m.b(i).el(k).fb;
            % Create the Lq vector for the beam 
            m.b(i).Lq(1,6*(k-1)+1:6*(k+1))=m.b(i).Lq(1,6*(k-1)+1:6*(k+1))+m.b(i).el(k).Lq;
            % Create Lb for the beam 
            m.b(i).Lb = m.b(i).Lb + m.b(i).el(k).Lb;    
            % Jx
            m.b(i).Jx = m.b(i).Jx + m.b(i).el(k).Jx; 
            % lp 
            m.b(i).lp = m.b(i).lp + m.b(i).el(k).lp; 
            % lb
            m.b(i).lb = m.b(i).lb + m.b(i).el(k).lb; 
            % lq
            m.b(i).lq(1,6*(k-1)+1:6*(k+1))=m.b(i).lq(1,6*(k-1)+1:6*(k+1))+m.b(i).el(k).lq;
            % Sq
            m.b(i).Sq(6*(k-1)+1:6*(k+1),1)=m.b(i).Sq(6*(k-1)+1:6*(k+1),1)+m.b(i).el(k).Sq;
            % fp
            m.b(i).fp(6*(k-1)+1:6*(k+1),1)=m.b(i).fp(6*(k-1)+1:6*(k+1),1)+m.b(i).el(k).fp;
        end
    end
end
m = m_compute_matrices(m);



