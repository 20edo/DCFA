
function el=el_init()
% This function initializes an element
% el.sc               Section of the element
% el.L                Length                        [mm]
% el.M                Mass matrix                   6*6 martix
% el.K                Stiffness matrix              6*6 matrix
% el.C                Dissipation matrix            6*6 matrix
% functions defined on el: 
% el_mass_assembly
% el_stiff_assembly
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
%% 
el.sc=sc_init();       % Section of the element
el.L=nan;       % Length                        [mm]
el.M=NaN(6);    % Mass matrix                   6*6 martix
el.K=NaN(6);    % Stiffness matrix              6*6 matrix
el.C=NaN(6);    % Dissipation matrix            6*6 matrix
el.fa = sparse(zeros(12,1));
el.Ka = sparse(zeros(12));
el.Ca = sparse(zeros(12));
el.fb = sparse(zeros(12,1)); 
el.Lq = sparse(zeros(1,12)); 
el.Lb = 0; 
el.Jx = 0; 
el.lp = 0; 
el.lb = 0; 
el.lq = sparse(zeros(1,12)); 
el.Sq = sparse(zeros(12,1));
el.fp = sparse(zeros(12,1));
el.Ham = sparse(zeros(12));
el.Ham_dk = sparse(zeros(12));
el.Hbam = sparse(zeros(12,1));
el.Hamb = sparse(zeros(1,12));
el.Hbb = 0; 
el.navier=sparse(zeros(12));    % Matrix used to recover the stresses in the correnti
el.load = sparse(zeros(12));
