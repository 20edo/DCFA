function m=m_add_load_beam(m,i,f)
% This function adds the load f to the i-th beam of the model
% m         model
% i         number of beam of reference
% f=@(x)    load                            [6*1]

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%

beam=m.b(i);
l=norm(m.b(i).in(2).x-m.b(i).in(1).x);      % Length of each element

% Definition of some useful functions
function y=f3(f,x)
    z=f(x);
    y=z(3);
end
function y=f2(f,x)
    z=f(x);
    y=z(3);
end

% First node
ftot=integral(@(x) f(x),beam.in(1).x,beam.in(1).x+l/2,'ArrayValued',true);                     % Total force on a node
mforce=integral(@(x) [0 -f3(f,x) f2(f,x)]'.*(beam.o+beam.vx*x),...
    beam.in(1).x,beam.in(1).x+l/2,'ArrayValued',true);                                         % Moment due to the force distribuition
mtot=mforce+ftot(4:6);                                                      % Total moment on a node

m.b(i).in(1).f=m.b(i).in(1).f+[ftot(1:3);mtot];                             % Add loads to the internal nodes
m.b(i).f(6*(1-1)+1:6*1)=m.b(i).f(6*(1-1)+1:6*1)+[ftot(1:3);mtot];           % Add loads to the beams


% General node
for j=2:(length(beam.in)-1)
    ftot=integral(@(x) f(x),beam.in(j).x-l/2,beam.in(j).x+l/2,'ArrayValued',true);             % Total force on a node
    mforce=integral(@(x) [0 -f3(f,x) f2(f,x)]'.*(beam.o+beam.vx*x),...
        beam.in(j).x-l/2,beam.in(j).x+l/2,'ArrayValued',true);                                 % Moment due to the force distribuition
    mtot=mforce+ftot(4:6);                                                  % Total moment on a node
    m.b(i).in(j).f=m.b(i).in(j).f+[ftot(1:3);mtot];                         % Add loads to the internal nodes
    m.b(i).f(6*(j-1)+1:6*j)=m.b(i).f(6*(j-1)+1:6*j)+[ftot(1:3);mtot];       % Add loads to the beams
end


% Last node
ftot=integral(@(x) f(x),beam.in(end).x-l/2,beam.in(end).x,'ArrayValued',true);                 % Total force on a node
mforce=integral(@(x) [0 -f3(f,x) f2(f,x)]'.*(beam.o+beam.vx*x),...
    beam.in(end).x-l/2,beam.in(end).x,'ArrayValued',true);                                     % Moment due to the force distribuition
mtot=mforce+ftot(4:6);                                                      % Total moment on a node
m.b(i).in(end).f=m.b(i).in(end).f+[ftot(1:3);mtot];                         % Add loads to the internal node
m.b(i).f(end-5:end)=m.b(i).f(end-5:end)+[ftot(1:3);mtot];                   % Add loads to the beam
    

end
