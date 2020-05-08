function m=m_add_beam(m,b)
% This function adds a beam to the model
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

% Create/assign first and last node to the beam
if any(m.en(:).x==b.o)
    node=en_free(b.on);
    m.en=[m.en node];
    b.on=length(m.en);
else
    b.on=find(m.en(:).x==b.o);
end

if any(m.en(:).x==b.o+b.L*b.v)
    node=en_free(b.o+b.L*b.v);
    m.en=[m.en node];
    b.en=length(m.en);
else
    b.en=find(m.en(:).x==b.o+b.L*b.v);
end

% Add beam to the model
m.b=[m.b b];
    