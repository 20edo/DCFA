function en=en_ground(position)
% This function creates a ground node
% Position---> 3*1 vector (position of the node)

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


if isequal(size(position),[1,3])
    position=position';
elseif ~isequal(size(position),[1,3])
    error('position is not defined properly');
end


% External node
en.x=position;              % Position in the 3d space [3*1]    [mm]
en.d=zeros(6,1);            % Displacement vector [6*1] vector  [mm]
en.c=true(6,1);              % Constraint 
en.M=zeros(6);                % Mass matrix
en.K=zeros(6);                % Stiffness matrix
en.C=zeros(6);                % Dissipation matrix
