function m=m_add_mass_forces(m,n,Z)
% This function adds to the model the mass forces.
% m-----> model
% n-----> load factor                   (default 1)
% Z-----> versor of the mass forces     (default [0 0 -1]')

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

%% Fix optional variables if not present
g=9.81;

if nargin==1
    n=1;
    Z=[0 0 -1]';
elseif nargin==2
    Z=[0 0 -1]';
end

for i=1:length(m.b)
    beam=m.b(i);
    l=norm(m.b(i).in(2).x-m.b(i).in(1).x);      % Length of each element
    vx=beam.vx;
    vy=beam.vy;
    vz=cross(vx,vy);
    vz=vz/norm(vz);

    % First node
    ftot=l/2*beam.el(1).sc.m*Z;                      % Total force on a node
    mtot=cross(l/2*beam.el(1).sc.m*Z,(l/4*vx+beam.el(1).sc.ycg*vy+beam.el(1).sc.zcg*vz));   % Moment due to cg!=(xnode,0,0)
%     mtot=[0 0 0]';
    m.b(i).in(1).f=m.b(i).in(1).f+g*n*[ftot(1:3);mtot];                       % Add loads to the internal nodes
    m.b(i).f(6*(1-1)+1:6*1)=m.b(i).f(6*(1-1)+1:6*1)+g*n*[ftot(1:3);mtot];     % Add loads to the beam


    % General node
    for j=2:(length(beam.in)-1)
        ftot=(l/2*beam.el(j).sc.m)*Z+l/2*beam.el(j-1).sc.m*Z;                      % Total force on a node
        mtot=cross(l/2*beam.el(j).sc.m*Z,(l/4*vx+beam.el(j).sc.ycg*vy+beam.el(j).sc.zcg*vz))+...      % "right" half-element and then "left " half-element
        cross(l/2*beam.el(j-1).sc.m*Z,(-l/4*vx+beam.el(j-1).sc.ycg*vy+beam.el(j-1).sc.zcg*vz));        % Moment due to cg!=(xnode,0,0)
%         mtot=[0 0 0]';
        m.b(i).in(j).f=m.b(i).in(j).f+g*n*[ftot(1:3);mtot];                           % Add loads to the internal node
        m.b(i).f(6*(j-1)+1:6*j)=m.b(i).f(6*(j-1)+1:6*j)+g*n*[ftot(1:3);mtot];         % Add loads to the beam
    end


    ftot=l/2*beam.el(end).sc.m*Z;                      % Total force on a node
    mtot= cross(l/2*beam.el(end).sc.m*Z,(-l/4*vx+beam.el(end).sc.ycg*vy+beam.el(end).sc.zcg*vz));      % "right" half-element and then "left " half-element % Moment due to cg!=(xnode,0,0)
%     mtot=[0 0 0]';
    m.b(i).in(end).f=m.b(i).in(end).f+g*n*[ftot(1:3);mtot];                           % Add loads to the internal node
    m.b(i).f(end-5:end)=m.b(i).f(end-5:end)+g*n*[ftot(1:3);mtot];                   % Add loads to the beam


end

end