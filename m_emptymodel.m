function m=m_emptymodel()
% This function generates an empty model
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
m.en=[];                    % Extrenal nodes (where more than one beam join or
                            % masses, forces, constraints,... are applied)
m.b=[];
m.M=nan;                    % Mass matrix of the model
m.K=nan;                    % Stiffness matrix of the model
