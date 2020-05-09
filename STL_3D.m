function [displ] = STL_3D(beam)

% This function gives a 3D representation of the displacements in
% u,v,w,theta.
% 
%
% Input:
% beam is the structure of the beam that contains all its proprierties

%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Pasturenzi Lorenzo    944610
%               Tacchi Alberto        944579
%               Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
%% Beam proprierties
L = beam.L;
nel = beam.nel;
%% Displacements
dL = beam.el(1).L;
displ = []; % vector of displacements [u v w th]'
for i = 1:nel-1
    xa = beam.in(i).x;
    xb = beam.in(i+1).x;
    clear x 
    x = linspace(xa, xb, 3); % we evaluate the displacements for 10 points inside each elements
    x = x(1:end-1); %otherwise we take into account 2 times the same node
    for h = 1:size(x,2)
        eps =(2*x(h)- (xb+xa))/(xb-xa);
        N = [(1 - eps)/2,                   0,                   0,         0,                             0,                            0, (1+eps)/2,                   0,                   0,         0,                              0,                             0;
             0, 1/4*(2-3*eps+eps^3),                   0,         0,                             0, 1/4*(1-eps-eps^2+eps^3)*dL/2,         0, 1/4*(2+3*eps-eps^3),                   0,         0,                              0, 1/4*(-1-eps+eps^2+eps^3)*dL/2;
             0,                   0, 1/4*(2-3*eps+eps^3),         0, -1/4*(1-eps-eps^2+eps^3)*dL/2,                            0,         0,                   0, 1/4*(2+3*eps-eps^3),         0, -1/4*(-1-eps+eps^2+eps^3)*dL/2,                             0;
             0,                   0,                   0, (1-eps)/2,                             0,                            0,         0,                   0,                   0, (1+eps)/2,                              0,                            0];
        displ_temp = N * [beam.in(i).d; beam.in(i+1).d];
        displ = [displ displ_temp];
        %displ = [u(x1) u(x2)...
        %         v(x1) v(x2)...
        %         w(x1) w(x2)...
        %         th(x1) th(x2)...]   4Xnel
    end
end

xa = beam.in(nel).x;
xb = beam.in(nel+1).x;
clear x 
x = linspace(xa, xb, 10); % now we are interested in the free end
displ_temp = N * [beam.in(i).d; beam.in(i+1).d];
displ = [displ displ_temp];

N_points = size(displ,2);
%% 3D plot
%% Create the beam structure
% These are geometric properties along [x,y,z]
W = L;
H = 1;
L = 1;
% Vertices of the beam
V = [0 0 0;
     W 0 0;
     W H 0;
     0 H 0;
     0 0 L;
     W 0 L;
     W H L;
     0 H L];
 % triangular facets of the beam (you need to start from somewhere)
 V = V - [0 H/2 L/2];
 F = [1 2 3;
      1 3 4;
      1 2 6;
      5 6 1;
      1 4 5;
      5 8 4;
      2 3 7;
      6 7 2;
      4 3 7;
      7 8 4;
      5 6 7;
      5 8 7];
patch('Faces',F, 'Vertices', V, 'FaceColor','blue', 'FaceAlpha',0.8 )
%% Save to STL file (not sure it is necessary...) 
% to later remesh it with a finer mesh
TR = triangulation(F,V);
stlwrite(TR,'wing.stl')
%% Import STL file and mesh it 3D
% we use tetraeda with only 4 vertexes (linear)
model = createpde(1);
importGeometry(model,'wing.stl');
OrigMesh = generateMesh(model, 'Hmin',L/N_points, 'Hmax', L/N_points, 'GeometricOrder', 'linear');
% if you don't close this figure the next one overlays on it
pdeplot3D(model)
%% Now we take into account the displacements
l = length(OrigMesh.Nodes);
PerturbNodes = OrigMesh.Nodes;
tmp = OrigMesh.Elements.';
u = displ(1,:);
v = displ(2,:);
w = displ(3,:);
ang = displ(4,:);

for j = 1:l
    Node = OrigMesh.Nodes(:,j);
    perturbNode = Node;
    k = round(Node(1)/(L/N_points));
    perturbNode(2) = Node(2)*cos(ang)-Node(3)*sin(ang);
    perturbNode(3) = Node(3)*cos(ang)+Node(2)*sin(ang); 
    perturbNode(1) = perturbNode(1)+u;
    perturbNode(2) = perturbNode(2)+v;
    perturbNode(3) = perturbNode(3)+w;
    PerturbNodes(:,j) = perturbNode;
end
    figure
    TR = triangulation(tmp,PerturbNodes.');
    [f,p] = freeBoundary(TR);
    trisurf(f,p(:,1),p(:,2),p(:,3), ...
       'FaceColor',rand(3,1),'FaceAlpha',0.8);
     axis('equal')
end




