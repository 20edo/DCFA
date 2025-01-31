function m=m_add_beam(m,b)
% This function adds a beam to the model
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

toll=1e-10;     % Distance between two node under which they are collapsed on the same node
% Create/assign first and last node to the beam
% Are the first or the node of the beam already in the model?
origin=false;
eend=false;
for i=1:length(m.en)
    if norm(m.en(i).x-b.o)<toll
        origin=i;
    elseif norm(m.en(i).x-(b.o+b.L*b.vx))<toll
        eend=i;
    end
end

if  origin       
    b.on=origin; % If there is a node coincident with the origin, save its coordinate
    
else            % else, create one
    node=en_free(b.o);
    m.en=[m.en node];
    b.on=length(m.en);
end

% Same for the end of the beam
if  eend
    b.en=eend;
else
    node=en_free(b.o+b.L*b.vx);
    m.en=[m.en node];
    b.en=length(m.en);
end

% Add beam to the model
m.b=[m.b b];
    