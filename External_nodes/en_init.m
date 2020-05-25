function en=en_init()
% This function initializes an external node
% 
% en.x            Position in the 3d space [3*1]    [mm]
% en.d            Displacement vector [6*1] vector  [mm]
% en.c            Constraint 
% en.f            External load                     [N]
% en.M            Mass matrix
% en.K            Stiffness matrix
% en.C            Dissipation matrix
% functions defined on en:
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
%%
en.x=nan(3,1);              % Position in the 3d space [3*1]    [mm]
en.d=zeros(6,1);            % Displacement vector [6*1] vector  [mm]
en.c=nan(6,1);              % Constraint 
en.f=zeros(6,1);            % External load                     [N]
en.M=nan(6);                % Mass matrix
en.K=nan(6);                % Stiffness matrix
en.C=nan(6);                % Dissipation matrix