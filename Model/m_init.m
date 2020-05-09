function m=m_init()
% This function initializes a model
% m.en                        % Extrenal nodes (where more than one beam join or
%                             % masses, forces, constraints,... are applied)
% m.b                         % Beams in the model
% m.M                         % Mass matrix of the model
% m.K                         % Stiffness matrix of the model
% 
% functions defined on model:
% m_add_beam


m.en=[];                    % Extrenal nodes (where more than one beam join or
                            % masses, forces, constraints,... are applied)
m.b=[];
m.M=nan;                    % Mass matrix of the model
m.C=0;                      % Dissipation matrix (defult value = 0)
m.K=nan;                    % Stiffness matrix of the model
