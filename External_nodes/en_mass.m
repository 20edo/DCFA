function en=en_mass(position,M)
% This function creates a lumped mass node
% Position---> 3*1 vector (position of the node)
% M       ---> mass or mass matrix

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


if isequal(size(position),[1,3])
    position=position';
elseif ~isequal(size(position),[3,1])
    error('position is not defined properly');
end

if isequal(size(M),[1 1])
    en.M=diag([M,M,M,0,0,0]);                % Mass matrix
elseif isequal(size(M),[6 6])
    en.M=M;                     % Mass matrix
else
    error('M is not defined properly');
end

% External node
en.x=position;              % Position in the 3d space [3*1]    [mm]
en.d=zeros(6,1);            % Displacement vector [6*1] vector  [mm]
en.c=false(6,1);              % Constraint 
en.f=zeros(6,1);            % External load                     [N]
en.K=zeros(6);                % Stiffness matrix
en.C=zeros(6);                % Dissipation matrix

