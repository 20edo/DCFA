function m_plot_eigenshape2(model,options,U)

% This function gives a 3D representation of the displacements in
% u,v,w,theta.
% 
%
% Input:
% beam is the structure of the beam that contains all its proprierties
% option is a structure including:
% plot = true / false - plot the beam with and without deformation
% plotColor = 'string' - face color of the plot
% saveSTL = true / false - save the beam to STL with and without deformation
% ovs = 4 - oversampling factor in the visualization of the displacements
% point_section - N/2 number of points in which i want to subdivide the half 
%                 wing section (2 points are enough for square
%                 sections)
% N = Number of modes
% alpha = shift for ROM_solver
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

%% Check
if (~exist('options.alpha', 'var'))
    % "dof" parameter does not exist
     options.alpha = 0.3;
end
if (~exist('options.saveImages', 'var'))
    % "dof" parameter does not exist
     options.saveImages = 0;
end
%% Parts of the model
for r = 1:size(model.b,2)
    beam = model.b(r);
    %% Beam proprierties
    dL = beam.el(1).L;
    L = beam.L;
    nel = beam.nel;
    %% 3D plot
    %% Create the beam structure
    % create the vector of the vertices & of the faces of the beam/wing/fuselage
    clear x
    V=[];
    F=[];
    
    if beam.cart == 1 %cartesian coordinates
        if beam.ssh
            n=options.point_section*2-2; % the 2 external points are the same for upper
            % section & bottom section
        else
        n = 4;
        end
        
        
        ymax = beam.Ymax(0);
        ymin = beam.Ymin(0);
        ymax_end = beam.Ymax(L);
        ymin_end = beam.Ymin(L);
        if beam.ssh
            angle_cgl=linspace(pi,pi/2,options.point_section);
            % Given the complex geometry of the wing we use the Chebyschev-Gaus_Lobatto
            % points to draw the section of the airfoil
            y_up=(cos(angle_cgl)+1)*(ymax-ymin);
            y_end_up=(cos(angle_cgl)+1)*(ymax_end-ymin_end);
        else
            y_up=[-0.5 0.5]*(ymax-ymin);
            y_end_up=[-0.5 0.5]*(ymax_end-ymin_end);
        end
        y_up=y_up';
        y_bot=flip(y_up);
        y_end_up=y_end_up';
        y_end_bot=flip(y_end_up);
        if beam.ssh
            y_bot=y_bot(2:end-1);
            y_end_bot=y_end_bot(2:end-1);
        end
        x_0_up=zeros(size(y_up,1),1);
        x_0_bot=zeros(size(y_bot,1),1);
        x_end_up=ones(size(y_end_up,1),1).*L;
        x_end_bot=ones(size(y_end_bot,1),1).*L;
        z_up=beam.Zmax;
        z_bot=beam.Zmin;
        if beam.ssh
            V=[x_0_up y_up z_up(x_0_up,y_end_up);
                x_0_bot y_bot z_bot(x_0_bot,y_end_bot);
                x_end_up y_end_up z_up(x_end_up,y_end_up);
                x_end_bot y_end_bot z_bot(x_end_bot,y_end_bot)];
        
            o = [zeros(size(y_up,1),1) beam.el(1).sc.yo*ones(size(y_up,1),1) zeros(size(y_up,1),1);
                zeros(size(y_bot,1),1) beam.el(1).sc.yo*ones(size(y_bot,1),1) zeros(size(y_bot,1),1);
                zeros(size(y_end_up,1),1) beam.el(end).sc.yo*ones(size(y_end_up,1),1) zeros(size(y_end_up,1),1);
                zeros(size(y_end_bot,1),1) beam.el(end).sc.yo*ones(size(y_end_bot,1),1) zeros(size(y_end_bot,1),1)];
            V= V + o;
        else
            V=[x_0_up y_up z_up(x_0_up,y_end_up);
                x_0_bot y_bot z_bot(x_0_bot,y_end_bot);
                x_end_up y_end_up z_up(x_end_up,y_end_up);
                x_end_bot y_end_bot z_bot(x_end_bot,y_end_bot)];
        end
        F = [1 2 n];
        for h=2:n/2 % faces on the section for x=0
            if h<n/2
                F=[F;
                    h h+1 n-h+2;
                    h+1 n-h+2 n-h+1];
            else %for h==n/2
                F=[F;
                    h h+1 n-h+2];
            end
        end
        for i=1:n-1 %faces on the upper & lower part of the span of the wing
            F=[F;
                i   i+1 i+n;
                i+1 i+n i+n+1];
        end
        F = [F;
            n n+1 2*n;
            1 n n+1];
        F=[F;
            n+1 n+2 2*n];
        for j=2:n/2 % faces on the section for x=L
            if j<n/2
                F=[F;
                    n+j n+j+1 2*n-j+2;
                    n+j+1 2*n-j+2 2*n-j+1];
            else %for j==n/2
                F=[F;
                    n+j n+j+1 2*n-j+2];
            end
        end
        
    elseif beam.pol == 1 % polar coordinates
        %        (4)o
        %         /   \
        %     (5)o     o(3)
        %      /         \
        %  (6)o     o(1)  o(2)
        %      \         /
        %     (7)o     o(9)
        %         \   /
        %        (8)o
        n=2*options.point_section;
        th_0=linspace(beam.Thmin(0),beam.Thmax(0),n+1);
        th_0=th_0(1:end-1);
        th_0=th_0';
        th_L=linspace(beam.Thmin(L),beam.Thmax(L),n+1);
        th_L=th_L(1:end-1);
        th_L=th_L';
        x_0=zeros(size(th_0,1),1);
        x_L=ones(size(th_L,1),1).*L;
        rho=beam.Rhomax;
        if rho(x_0,th_0) == 0
            x_0 = [];
            y_0 = [];
            z_0 = [];
        else
            y_0=rho(x_0,th_0).*cos(th_0);
            z_0=rho(x_0,th_0).*sin(th_0);
        end
        y_L=rho(x_L,th_L).*cos(th_L);
        z_L=rho(x_L,th_L).*sin(th_L);
        V=[0 0 0;
            x_0 y_0 z_0;
            L 0 0;
            x_L y_L z_L];
        if rho(0,th_0) == 0
            for i=1:n-1 %faces on the lateral area of the cylinder
                F=[F;
                    1   i+2 i+3];
            end
            F = [F;
                1 3 n+2];
            
            for j=1:n % faces on the section for x=L
                if j<n
                    F=[F;
                        2 j+2 j+3];
                else %for j==n
                    F=[F;
                        2 n+2 3];
                end
            end
        else
            for h=1:n % faces on the section for x=0
                if h<n
                    F=[F;
                        1 h+1 h+2];
                else %for h==n
                    F=[F;
                        1 n+1 2];
                end
            end
            for i=1:n-1 %faces on the lateral area of the cylinder
                F=[F;
                    i+1   i+2 i+n+2;
                    i+2 i+n+2 i+n+3];
            end
            F = [F;
                n+1 2 2*n+2;
                2 2*n+2 n+3];
            
            for j=1:n % faces on the section for x=L
                if j<n
                    F=[F;
                        n+2 n+j+2 n+j+3];
                else %for j==n
                    F=[F;
                        n+2 2*n+2 n+3];
                end
            end
        end
        
    end
    vx_l = beam.vx; % axis versor of the beam (local reference)
    vy_l = beam.vy;
    vz_l = cross(vx_l,vy_l);
    vx_g = [1,0,0]';  % axis versor of the beam (global reference)
    vy_g = [0,1,0]';
    vz_g = [0,0,1]';
    T = zeros(3);
    T(1,1) = dot(vx_g, vx_l);
    T(1,2) = dot(vx_g, vy_l);
    T(1,3) = dot(vx_g, vz_l);
    T(2,1) = dot(vy_g, vx_l);
    T(2,2) = dot(vy_g, vy_l);
    T(2,3) = dot(vy_g, vz_l);
    T(3,1) = dot(vz_g, vx_l);
    T(3,2) = dot(vz_g, vy_l);
    T(3,3) = dot(vz_g, vz_l);
    V = T*V';
    V = V';
    %% Save the mesh in the struct Support
%     TR = triangulation(F,V);
%     fname = sprintf('%s', beam.name);
%     stlwrite(TR,'fname')
    modello = createpde(1);
    geometryFromMesh(modello,V',F');
    % Mesh the beam with 'linear' tetrahedra, to save memory and complexity.
    OrigMesh = generateMesh(modello, ...
        'Hmin',L/nel,'Hmax', L/nel, ... 
        'GeometricOrder', 'linear');
    Nodes = OrigMesh.Nodes + beam.o;
    tmp = OrigMesh.Elements.';
    Support(r).Initial_Mesh= OrigMesh;
    Support(r).Mesh_elements= tmp;
    Support(r).Initial_nodes= Nodes;
end

clear r h i j
%% I want to find the firsts N modes


%% Add the zeros for the constrained displacements
for j=1:length(model.en)
    for k=1:6
        if model.en(j).c(k) 
            index = 6*(j-1)+k;
            U = [U(1:index-1,:);zeros(1,size(U,2));U(index:end,:)]; 
        end
    end            
end
%% Add the eigenshapes into the displacements
for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*U; 
    for k = 1:model.b(i).nel+1
       model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end
clear i k
for r = 1:size(model.b,2)
    beam = model.b(r);
    %% Beam proprierties
    dL = beam.el(1).L;
    L = beam.L;
    nel = beam.nel;
    %% Beam rotation
    vx_l = beam.vx; % axis versor of the beam (local reference)
    vy_l = beam.vy;
    vz_l = cross(vx_l,vy_l);
    vx_g = [1,0,0]';  % axis versor of the beam (global reference)
    vy_g = [0,1,0]';
    vz_g = [0,0,1]';
    T = zeros(3);
    T(1,1) = dot(vx_g, vx_l);
    T(1,2) = dot(vx_g, vy_l);
    T(1,3) = dot(vx_g, vz_l);
    T(2,1) = dot(vy_g, vx_l);
    T(2,2) = dot(vy_g, vy_l);
    T(2,3) = dot(vy_g, vz_l);
    T(3,1) = dot(vz_g, vx_l);
    T(3,2) = dot(vz_g, vy_l);
    T(3,3) = dot(vz_g, vz_l);
    % assembly of the big rotation matrix
    N = size(beam.M,1)/3;
    Ar = repmat(T, 1, N);
    Ac = mat2cell(Ar, size(T,1), repmat(size(T,2),1,N));
    R = sparse(blkdiag(Ac{:}));
    clear N
    %% Now we take into account the displacements
    clear i
    clear j
    OrigMesh = Support(r).Initial_Mesh;
    l = length(OrigMesh.Nodes);
    PerturbNodes = OrigMesh.Nodes;
    for n = 1:options.N
        for j = 1:l
            Node                = OrigMesh.Nodes(:,j);
            perturbNode         = Node;
            % Nodes in local frame
            LocalNode = T'*perturbNode;
            x                   = LocalNode(1);
            el                  = (floor(x/dL)+1);
            if el > nel %for points on the last section
                el = nel;
            end
            if el <= 0 % this is because there are elements with x=-2.2138e-16
                el = 1;
            end
            xa                  = beam.in(el).x;
            xb                  = beam.in(el+1).x;
            eps                 = (2*x- (xb+xa))/(xb-xa);
            N                   = [(1 - eps)/2,                   0,                   0,         0,                             0,                            0, (1+eps)/2,                   0,                   0,         0,                              0,                             0;
                0, 1/4*(2-3*eps+eps^3),                   0,                   0,                             0, 1/4*(1-eps-eps^2+eps^3)*dL/2,         0, 1/4*(2+3*eps-eps^3),                   0,         0,                              0, 1/4*(-1-eps+eps^2+eps^3)*dL/2;
                0,                   0, 1/4*(2-3*eps+eps^3),                   0, -1/4*(1-eps-eps^2+eps^3)*dL/2,                            0,         0,                   0, 1/4*(2+3*eps-eps^3),         0, -1/4*(-1-eps+eps^2+eps^3)*dL/2,                             0;
                0,                   0,                   0,           (1-eps)/2,                             0,                            0,         0,                   0,                   0, (1+eps)/2,                              0,                             0;
                0,                   0, -1/4*(-3+3*eps^2)*2/dL,                   0,            1/4*(-1-2*eps+3*eps^2),                            0,         0,                   0, -1/4*(+3-3*eps^2)*2/dL,         0,             1/4*(-1+2*eps+3*eps^2),                             0;
                0,  1/4*(-3+3*eps^2)*2/dL,                   0,                   0,                             0,           1/4*(-1-2*eps+3*eps^2),         0,  1/4*(+3-3*eps^2)*2/dL,                   0,         0,                              0,           1/4*(-1+2*eps+3*eps^2)];
            d = R(1:12,1:12)'*[beam.in(el).d(:,n); beam.in(el+1).d(:,n)];% displacements in local frame
            displ               = N * d;
            % Initialize local perturbation variable
            LocalPerturbNode=LocalNode;
            % Rotating the section in the local frame
            LocalPerturbNode(2)      = LocalNode(2)*cos(displ(4))-LocalNode(3)*sin(displ(4));
            LocalPerturbNode(3)      = LocalNode(3)*cos(displ(4))+LocalNode(2)*sin(displ(4));
            % Applying the displacements in the local frame
            LocalPerturbNode(1)      = LocalPerturbNode(1)+displ(1);
            LocalPerturbNode(2)      = LocalPerturbNode(2)+displ(2);
            LocalPerturbNode(3)      = LocalPerturbNode(3)+displ(3);
            % Back to the global reference
            perturbNode = T*LocalPerturbNode;
            PerturbNodes(:,j)   = perturbNode;
            
        end
        PerturbNodes = PerturbNodes + beam.o;
        Support(r).Deformed_nodes(n).modes = PerturbNodes;
    end
end

%% Plot deformed configuration
if options.plot_deformed
    for i = 1:options.N
        for r = 1:size(model.b,2)
            TR = triangulation(Support(r).Mesh_elements,Support(r).Deformed_nodes(i).modes.');
            [f,p] = freeBoundary(TR);
            trisurf(f,p(:,1),p(:,2),p(:,3), ...
                'FaceColor',options.plotColor,'FaceAlpha',0.8);
            axis('equal');
%             title({['Deformed model for mode:',num2str(i)];['Frequency:',num2str(w)]});
            hold on
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        if options.saveImages
            for h = 1:3
                f = figure('visible','off');
                if h == 1
                    view([1 0 0])
                    fname = sprintf('Mode\t%d\tx', i);
                    saveas(f,'fname','svg')
                elseif h==2
                    view([0 1 0])
                    fname = sprintf('Mode\t%d\ty', i);
                    saveas(f,'fname','svg')
                elseif h==3
                    view([0 0 1])
                    fname = sprintf('Mode\t%d\tz', i);
                    saveas(f,'fname','svg')
                end
            end
        end
    end
end

if options.saveSTL
    TR = triangulation(tmp,PerturbNodes.');
    [f,p] = freeBoundary(TR);
    TR = triangulation(f,p);
    stlwrite(TR,'beamDeformed.stl')
end

end