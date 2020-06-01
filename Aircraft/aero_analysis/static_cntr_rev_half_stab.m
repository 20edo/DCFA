% Static aero analysis of the clamped wing
% - Control reversal
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
clear all, close all, clc
cd ..
%% Generate the wing model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\


%% Build the swept wing model
half_stab=m_init();
half_stab.en=[en_ground(aircraft.en(14).x)];
stab_list=[aircraft.b(15)];
for i=1:length(stab_list)
    half_stab=m_add_beam(half_stab,stab_list(i));
end
%% Generate the straight wing model
wing_straight = half_stab;
%% Add the aero loads
half_stab = m_add_aero_loads(half_stab,[1,0,0]');
wing_straight = m_add_aero_loads_straight(wing_straight,[1,0,0]');


%% Assebly stiffness the K matrices for cntr_rev
K_cr = [half_stab.K, zeros(size(half_stab.K,1),1); zeros(1,size(half_stab.K,2)),0]; 
Ka_cr = [half_stab.Ka, half_stab.fb; half_stab.Lq, half_stab.Lb]; 

%% Find the control reversal (cr) dynamic pressure
[V_cr,D_cr]= eig(full(K_cr),full(Ka_cr));
q_cr = diag(D_cr);
[q_cr,I] = sort(real(q_cr));
V_cr = V_cr(:,I);                         % sort the eigenshapes
V_cr(:,q_cr<1)=[];                    % select the eigenshapes with positive eig
q_cr(q_cr<1)=[];
q_cr = q_cr(1);                   % select the minimum positive q_inf
V_cr = V_cr(:,1);                     % select its eigenshape

[V,D]= eig(full(half_stab.K),full(half_stab.Ka));
q_div = diag(D);
[q_div,I] = sort(real(q_div));
V = V(:,I);                             % sort the eigenshapes
V(:,q_div<0)=[];                        % select the eigenshapes with positive eig
q_div(q_div<0)=[];
q_div = q_div(1);                       % select the minimum positive q_inf
V_div = V(:,1);  

    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    m_plot_eigenshape(half_stab,options,V_div/20)
    
    