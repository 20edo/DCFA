function in=in_init()
% Initialize an internal node
% in.x                    Position of the node along the beam axis    [mm]
% in.d                    Displacement of the node [6*1]              [mm]
% in.f                    External load                               [N]

% functions defined on internal node;
%
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
in.x=nan;
in.d=zeros(6,1);
in.f=zeros(6,1);            % External load                     [N]
