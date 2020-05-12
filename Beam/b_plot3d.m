function [fig] = b_plot3d(beam,options)

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
    x = linspace(xa, xb, options.ovs); % we evaluate the displacements for 10 points inside each elements
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
        %         th(x1) th(x2)...]   4XN_points
    end
end

xa = beam.in(nel).x;
xb = beam.in(nel+1).x;
clear x 
x = linspace(xa, xb, options.ovs); % now we are interested in the free end
for h = 1:size(x,2)
    eps =(2*x(h)- (xb+xa))/(xb-xa);
    N = [(1 - eps)/2,                   0,                   0,         0,                             0,                            0, (1+eps)/2,                   0,                   0,         0,                              0,                             0;
         0, 1/4*(2-3*eps+eps^3),                   0,         0,                             0, 1/4*(1-eps-eps^2+eps^3)*dL/2,         0, 1/4*(2+3*eps-eps^3),                   0,         0,                              0, 1/4*(-1-eps+eps^2+eps^3)*dL/2;
         0,                   0, 1/4*(2-3*eps+eps^3),         0, -1/4*(1-eps-eps^2+eps^3)*dL/2,                            0,         0,                   0, 1/4*(2+3*eps-eps^3),         0, -1/4*(-1-eps+eps^2+eps^3)*dL/2,                             0;
         0,                   0,                   0, (1-eps)/2,                             0,                            0,         0,                   0,                   0, (1+eps)/2,                              0,                            0];
    displ_temp = N * [beam.in(nel).d; beam.in(nel+1).d];
    displ = [displ displ_temp];
end

N_points = size(displ,2);
%% 3D plot
%% Create the beam structure
% Create the beam structure.
% These are geometric properties along [x,y,z]
% Assumes the section of the beam to be 1/10 of its length
W = L;
H = W/10;
LL = W/10;
% Vertices of the beam
V = [0 0 0;
     W 0 0;
     W H 0;
     0 H 0;
     0 0 LL;
     W 0 LL;
     W H LL;
     0 H LL];
 % Center beam around (0,0)
 V = V - [0 H/2 LL/2];
 % Triangular facets of the beam (you need to start from somewhere)
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
% Save to STL file
if options.saveSTL
    TR = triangulation(F,V);
    stlwrite(TR,'beam.stl')
end

%% Mesh the beam in 3D
% we use tetraheda with only 4 vertexes (linear)
model = createpde(1);
geometryFromMesh(model,V.',F.');
% Mesh the beam with 'linear' tetrahedra, to save memory and complexity.
% In any case the mesh is not used to compute the solution
% The element size is the same for which the displacements have been computed.
OrigMesh = generateMesh(model, ...
    'Hmin',L/(N_points-1), 'Hmax', L/(N_points-1), ...
    'GeometricOrder', 'linear');
if options.plot
    TR = triangulation(OrigMesh.Elements.',OrigMesh.Nodes.');
    [f,p] = freeBoundary(TR);
    trisurf(f,p(:,1),p(:,2),p(:,3), ...
        'FaceColor',options.plotColor,'FaceAlpha',0.8);
    %pdeplot3D(model)
    axis('equal');
    title('Original beam without deformation');
end
%% Now we take into account the displacements
l = length(OrigMesh.Nodes);
PerturbNodes = OrigMesh.Nodes;
tmp = OrigMesh.Elements.';
u = displ(1,:);
v = displ(2,:);
w = displ(3,:);
ang = displ(4,:);

for j = 1:l
    Node                = OrigMesh.Nodes(:,j);
    perturbNode         = Node;
    k                   = round(Node(1)/(L/(N_points-1))) + 1;
    perturbNode(2)      = Node(2)*cos(ang(k))-Node(3)*sin(ang(k));
    perturbNode(3)      = Node(3)*cos(ang(k))+Node(2)*sin(ang(k));
    perturbNode(1)      = perturbNode(1)+u(k);
    perturbNode(2)      = perturbNode(2)+v(k);
    perturbNode(3)      = perturbNode(3)+w(k);
    PerturbNodes(:,j)   = perturbNode;
end
if options.plot
    fig = figure
    TR = triangulation(tmp,PerturbNodes.');
    [f,p] = freeBoundary(TR);
    trisurf(f,p(:,1),p(:,2),p(:,3), ...
        'FaceColor',options.plotColor,'FaceAlpha',0.8);
    axis('equal');
    title('Deformed beam');
end
if options.saveSTL
    TR = triangulation(tmp,PerturbNodes.');
    [f,p] = freeBoundary(TR);
    TR = triangulation(f,p);
    stlwrite(TR,'beamDeformed.stl')
end
end



