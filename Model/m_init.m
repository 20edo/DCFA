function m=m_init()
% This function initializes a model
% m.en                        % Extrenal nodes (where more than one beam join or
%                             % masses, forces, constraints,... are applied)
% m.b                         % Beams in the model
% m.M                         % Mass matrix of the model
% m.K                         % Stiffness matrix of the model
% m.f                         % External loads vector
% 


% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%

m.en=[];                    % Extrenal nodes (where more than one beam join or
                            % masses, forces, constraints,... are applied)
m.b=[];
m.M=nan;                    % Mass matrix of the model
m.C=0;                      % Dissipation matrix (defult value = 0)
m.f=[];                     % External loads vector

m.K=nan;                    % Stiffness matrix of the model
