%% *Wing Mounted engine (DCFA D.4.3)* 
%% Data initialization 

clear; clc; close all
%% 
% Wing span [m]

b = 20;
%% 
% Bending stiffness [N*m^2]

EJ = 1e6; 
%% 
% Torsional stiffness [N*m^2]

GJ = 1e4; 
%% 
% Mass per unit span [Kg/m]

m = 50; 
%% 
% Polar inertia per unit span, [Kg*m^2/m]

J_p = 10; 
%% 
% Center of mass offset, [m]

d = 0.2; 
%% 
% Mass of engine [kg]

M_engine = 500; 
%% 
% Forward offset of engine CM, [m]

y_engine = 1; 
%% 
% Vertical offset of engine CM, [m]

z_engine = -.5; 
%% 
% Inertia tensor of engine, [Kg*m^2]

J_engine = diag([50, 1, 50]); 
%% 
% Number of elements on the wing

N = 20; 
%% 
% Engine's node

N_engine = 8
%% 
% Element length

DL = b / N;
%% 
% Change number of modes

i_ROM_max =9
%% Build mass matrix and stiffness matrix
% $$0 = \delta \mathcal{W} = -\int_l \delta w_{/xx}^T EJ w_{/xx} dx - \int_l 
% \delta\theta_{/x}^TGJ\theta_{/x}dx-\int_l\delta_{CM}^Tm\ddot_{w}_{CM}dx - \int_{l}\delta\theta^TJ_P\ddot{\theta} 
% dx + \delta \mathcal{W}_f$$
% 
% $$O = \delta\mathcal{W} = \delta \mathbf{u}^T(-\mathbf{M}\mathbf{\ddot{u}} 
% - \mathbf{Ku}+\mathbf{f}(t))$$

Meww = m * DL / 420 * [156, -22 * DL, 54, 13 * DL;
                        -22 * DL, 4 * DL^2, -13 * DL, -3 * DL^2;
                        54, -13 * DL, 156, 22 * DL;
                        13 * DL, -3 * DL^2, 22 * DL, 4 * DL^2]; 
Keww = EJ / DL^3 * [12, -6 * DL, -12, -6 * DL;
                    -6 * DL, 4 * DL^2, 6 * DL, 2 * DL^2;
                    -12, 6 * DL, 12, 6 * DL;
                    -6 * DL, 2 * DL^2, 6 * DL, 4 * DL^2];
Mett = (J_p + m * d^2) * DL * [1/3, 1/6;
                                1/6, 1/3];
Kett = GJ / DL * [1, -1;
                -1, 1];
Mewt = -d * m * DL / 60 * [21, 9;
                        -3 * DL, -2 * DL;
                        9, 21;
                        2 * DL, 3 * DL];
%% Organize degrees of freedom as {$w_i$ , $\varphi_i$, $\theta_i$}

Me([3, 6], [3, 6]) = Mett;
Me([1, 2, 4, 5], [1, 2, 4, 5]) = Meww;
Me([1, 2, 4, 5], [3, 6]) = Mewt;
Me([3, 6], [1, 2, 4, 5]) = Mewt';

Ke([3, 6], [3, 6]) = Kett;
Ke([1, 2, 4, 5], [1, 2, 4, 5]) = Keww;

%% 
% Whole problem matrices

M = zeros(3 * (N + 1), 3 * (N + 1));
K = zeros(3 * (N + 1), 3 * (N + 1));

%% 
% Assembly

for i = 1:N
    M(3 * (i-1) +1:3 * (i+1),3 * (i-1) +1:3 * (i+1)) = M(3 * (i-1) +1:3 * (i+1),3 * (i-1) +1:3 * (i+1)) + Me;
    K(3 * (i-1) +1:3 * (i+1),3 * (i-1) +1:3 * (i+1)) = K(3 * (i-1) +1:3 * (i+1),3 * (i-1) +1:3 * (i+1)) + Ke;
end
% assembling matrix 6x6 considering 3X3 element in common
%% 
% Clamp constraint at the free end

M = M(4:end, 4:end); % for the free wing
K = K(4:end, 4:end);
%% 
% Engine contribution: kinematics 


Z_engine = [0, z_engine, 0;
            0 , 0, -z_engine;
            1, 0, y_engine]; % see [D.232]
% Optional: Add the inertia term of the motor to the mass matrix
MM_engine = Z_engine' * M_engine * Z_engine;

%% 
% Eigenvalues of the clamped wing 


[V, w2] = eig(K, M);
w2 = diag(w2);
[~, II] = sort(w2); % we sort to the lower to the higher
w = sqrt(w2(II));
V = [zeros(3, 3 * N); V(:, II)];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);
end


%% 
% Node coordinates (b is the wing span) 

x = [0:N]' / N * b;
%% 
% Chose how many modes we want to keep in the ROM


%i_ROM = [1:15]';

i_ROM = [1:i_ROM_max]';

%%  Calculation of the "exact" modes with engine 
% $$(\mathbf{M}+\mathbf{M}_{\mathrm{engine}})\mathbf{\ddot{u}}+\mathbf{Ku} = 
% \mathbf{0}$$
% 
% Insert the engine contribution from in the corresponding dof of the mass matrix. 

MM = M; %we are adding to the mass matrix of the "clear wing" the mass matrix of the engine
MM(3 * N_engine - 2:3 * N_engine, 3 * N_engine - 2:3 * N_engine) = MM(3 * N_engine - 2:3 * N_engine, 3 * N_engine - 2:3 * N_engine) + MM_engine; 
%% 
% Eig of the exact solution. 


[V, w2] = eig(K, MM);
w2 = diag(w2);
[dmy, II] = sort(w2);
w = sqrt(w2(II));
V = [zeros(3, 3 * N); V(:, II)];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i); % normalization 
end

V_exact = V(:, i_ROM);
w_exact = w(i_ROM); % vector of the eigenmode i wanted to consider
% reduce the matrix
K_exact = V_exact(4:end, :)' * K * V_exact(4:end, :); % I'm subtracting the eigenshapes related to the constrained nodes (see constrained mass matrix & stiffness matrix)
M_exact = V_exact(4:end, :)' * MM * V_exact(4:end, :);
%% Corrected a posteriori
% Calculate the exact modes of the clean wing, then insert the engine contribution 
%% 
% # $\mathbf{M}\mathbf{U}_{cl}\mathrm{Diag}\{\omega_{\mathrm{cl}}^2\}= \mathbf{KU}_{\mathrm{cl}}$
% # $\mathbf{s}_\mathrm{engine} \approx \mathbf{ZU}_{\mathrm{cl}}\mathbf{q}$
% # $(\mathrm{Diag\{m_\mathrm{cl}\}+\mathbf{m}_{\mathrm{engine}}} )\mathbf{\ddot{q}} 
% + \mathrm{Diag}\{k_{cl}\}\mathbf{q} = \mathbf{U}_{cl}^T\mathbf{f} $

[V, w2] = eig(K, M);
w2 = diag(w2);
[dmy, II] = sort(w2);
w = sqrt(w2(II));
V = [zeros(3, 3 * N); V(:, II)];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2),
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i); % normalization with max abs. 
end

V_corrected = V(:, i_ROM); % take the reduced base 
%insert the engine contribution and on the reduced matrix 
K_corrected = V_corrected(4:end, :)' * K * V_corrected(4:end, :); % (4:end) to take out the clamp constrain as before
M_corrected = V_corrected(4:end, :)' * M * V_corrected(4:end, :) + ...
                V_corrected(3 * N_engine + 1:3 * N_engine + 3, :)' * MM_engine * V_corrected(3 * N_engine + 1:3 * N_engine + 3, :);
                % contibution of the engine in the "clear" system using teh "clear" modeshapes

[VV, ww2] = eig(K_corrected, M_corrected);
ww2 = diag(ww2);
[dmy, II] = sort(ww2);
w_corrected = sqrt(ww2(II));
V = V_corrected * VV(:, II);
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2),
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);

    if (V(:, i)' * V_exact(:, i) < 0),
        V(:, i) = -V(:, i);
    end

end

V_corrected = V;
%% Static shapes
% First solve the static problem: $\mathbf{Ku} = \mathbf{Sf}_p$

[V, w2] = eig(K, M);
w2 = diag(w2);
[dmy, II] = sort(w2);
w = sqrt(w2(II));
V = V(:, II);
% solve for the static problem
U = K\[zeros(3 * N_engine - 3, 3); eye(3); zeros(size(K, 1) - 3 * N_engine, 3)];
% All zero matrix apart from the identity matrix corresponding to the
% location in which there is the engine
%% 
% $$\mathbf{U}_{\mathrm{ROM}} = [\mathbf{U}_{cl} \ \mathbf{U}_{st}]$$

V = [V(:, i_ROM(1:end - 3)), U]; % i add the static solution at the shapemode of the clear system
V = [zeros(3, size(V, 2)); V];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];
for i = 1:size(V, 2),
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i)./V(ij, i);
end

V_static = V;

%% 
% $$\mathbf{K}_{\mathbf{static}} = \mathbf{U_\mathrm{ROM}}^T\mathbf{K}\mathbf{U}_{\mathrm{ROM}}$$

K_static = V_static(4:end, :)' * K * V_static(4:end, :); 

%% 
% $$\mathbf{M}_{\mathrm{static}} = \mathbf{U}^T_{\mathrm{ROM}}(\mathbf{M} + 
% \mathbf{M}_{\mathrm{engine}}) \mathbf{U}_{\mathrm{ROM}} $$

M_static = V_static(4:end, :)' * M * V_static(4:end, :) +...
            V_static(3 * N_engine + 1:3 * N_engine + 3, :)' * MM_engine * V_static(3 * N_engine + 1:3 * N_engine + 3, :);
            % as did before (fot the corrected "a posteriori" method)
[VV, ww2] = eig(K_static, M_static);
ww2 = diag(ww2);
[dmy, II] = sort(ww2);
w_static = sqrt(ww2(II));
V = V_static*VV(:, II);
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];
for i = 1:size(V, 2),
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i)./V(ij, i);
    if (V(:, i)'*V_exact(:, i) < 0),
        V(:, i) = -V(:, i);
    end
end
V_static = V;

%% Substructuring: Craig and Bamptom 
% Partition the dof of the problem in reduced (r) and interface (i). 

i_i = [3 * N_engine - 2:3 * N_engine]'; % interface dof (in contact with the engine)
i_r = [[1:3 * N_engine - 3], [3 * N_engine + 1:3 * N]]'; % reduced dof (not in contact with the engine)
KK = K(i_r, i_r); 
MM = M(i_r, i_r);
%% 
% $$\mathbf{M}_{rr}\ddot{\mathbf{u}}_r +  \mathbf{K}_{rr}\ddot{\mathbf{u}}_r  
% = 0$$

[VV, ww2] = eig(KK, MM);
V = [VV([1:3*N_engine - 3], :); zeros(3, 3 * N - 3); VV(3 * N_engine - 2:end, :)];
%% 
% Consider the reduced equation: $\mathbf{K_{rr}u_{r}}+\mathbf{K_{ri}u_i} = 
% 0$, we have that $\mathbf{u_r} = -\mathbf{K_{rr}}^{-1}\mathbf{K_{ri}u_i}$

% we are now solving the equivalent of the static reduced problem (f=0 & u''=0)
U = -KK \ K(i_r, i_i); 

U = [U(1:3 * N_engine - 3, :); eye(3); U(3 * N_engine - 2:end, :)]; 
% as did before for the static shapes: we are adding an identity matrix in corresponding of the
% nodes of the engine's dof

%% 
% $U_{\textrm{cb}}$ = $\left\lbrack \begin{array}{cc}U_r  & -{{K_{\left\lbrack 
% \textrm{rr}\right\rbrack } }^{-1} K}_{\left\lbrace \textrm{ri}\right\rbrace 
% } \\0 & I\end{array}\right\rbrack$

V = [V(:, i_ROM(1:end -3)), U]; % we add the new solution at the shapemode of the clear system
V = [zeros(3, size(V, 2)); V];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);
end

V_cb = V;  
K_cb = V_cb(4:end, :)' * K * V_cb(4:end, :); 
M_cb = V_cb(4:end, :)' * M * V_cb(4:end, :) +...
        V_cb(3 * N_engine + 1:3 * N_engine + 3, :)' * MM_engine * V_cb(3 * N_engine + 1:3 * N_engine + 3, :);

[VV, ww2] = eig(K_cb, M_cb);
ww2 = diag(ww2);
[dmy, II] = sort(ww2);
w_cb = sqrt(ww2(II));
V = V_cb * VV(:, II);
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]']; 

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);

    if (V(:, i)' * V_exact(:, i) < 0)
        V(:, i) = -V(:, i);
    end

end

V_cb = V;

%%  Boundary mass modes with engine
% Insert a very large mass to simulate a clamped constrain in correspondance 
% of the dof of the engine 

MM = M; 
MM(3 * N_engine - 2:3 * N_engine, 3 * N_engine - 2:3 * N_engine) = MM(3 * N_engine - 2:3 * N_engine, 3 * N_engine - 2:3 * N_engine) + 1e6 * eye(3);
% we are adding a very large mass at the nodes (dof) corresponding to the engine
% == BOUNDARY MASS

% we used 1e6, but three order of magnitute above to the mass of the clear
% wing is enough to obtain a very tiny error

%% 
% $$(\mathbf{M}+\mathbf{M}_{bm})\mathbf{U}_{bm} \mathrm{Diag}\{\omega_{bm}^2\} 
% = \mathbf{K}\mathbf{U}_{bm}$$

[V, w2] = eig(K, MM);
w2 = diag(w2);
[dmy, II] = sort(w2);
w = sqrt(w2(II));
V = [zeros(3, 3 * N); V(:, II)];
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);
end

V_bm = V(:, i_ROM);
% insert the true engine mass matrix in the ROM model
% now that we have solved the eigenproblem, we use the REAL mass matrix &
% mass engine to find the solutions
%(we use the eigenshape of the fictitiuos mass to obtain results with very
%little error)
K_bm = V_bm(4:end, :)' * K * V_bm(4:end, :);
M_bm = V_bm(4:end, :)' * M * V_bm(4:end, :) +...
        V_bm(3 * N_engine + 1:3 * N_engine + 3, :)' * MM_engine * V_bm(3 * N_engine + 1:3 * N_engine + 3, :);
[VV, ww2] = eig(K_bm, M_bm);
ww2 = diag(ww2);
[dmy, II] = sort(ww2);
w_bm = sqrt(ww2(II));
V = V_bm * VV(:, II);
ii = [[1:3:3*(N+1)]'; [3:3:3*(N+1)]'];

for i = 1:size(V, 2)
    ij = find(abs(V(:, i)) == max(abs(V(ii, i))));
    V(:, i) = V(:, i) ./ V(ij, i);

    if (V(:, i)' * V_exact(:, i) < 0)
        V(:, i) = -V(:, i);
    end

end

V_bm = V;

% compare results...
disp('mode        ''exact''         corrected       static          Craig-Bampton   boundary mass');

for i = 1:length(w_exact)
    fprintf('%4d    %16f%16f%16f%16f%16f\n', i, w_exact(i), w_corrected(i), w_static(i), w_cb(i), w_bm(i));
end

for i = 1:length(i_ROM)
    figure;
    subplot(2, 1, 1);
    plot(x, V_exact(1:3:end, i), x, V_corrected(1:3:end, i), x, V_static(1:3:end, i), x, V_cb(1:3:end, i), x, V_bm(1:3:end, i));
    ylabel('w');
    title(sprintf('mode %d: \\omega=%e', i, w_exact(i)));
    subplot(2, 1, 2);
    plot(x, V_exact(3:3:end, i), x, V_corrected(3:3:end, i), x, V_static(3:3:end, i), x, V_cb(3:3:end, i), x, V_bm(3:3:end, i));
    xlabel('x');
    ylabel('\theta');
    legend('exact', 'corrected', 'static', 'Craig-Bampton', 'boundary mass',"Location","bestoutside");
    
end