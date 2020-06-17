% Convergence study of the clamped wing

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

%% Generate the aircraft model
cd ..
cd ..
cd generate_model\
generate_model

%% Load the exact solution

cd convergence_studies\Clamped_wing\ % Move to the right folder 

load w_esatta.mat
load V_esatti.mat

cd ..
cd ..
cd generate_model
%% Perform analysis

nel_tot=30:15:300;
t=zeros(size(nel_tot));
number=30;      % Number of eigenvalues considered
w_esatta(number+1:end)=[];
V_esatti(:,number+1:end)=[];

for i =1:length(nel_tot)
    tic
    % Build model
    model=build_clamped_wing(aircraft,nel_tot(i));
    model=m_compute_matrices(model);
    % Shift
    alpha=0;
    % Solve
    [V,D,flag] = eigs(model.K+alpha*model.M,model.M,number,'smallestabs');
    w=diag(D-alpha).^0.5;
    t(i)=toc;
    % Calculate errors
    error1=norm(w-w_esatta,1)/norm(w_esatta,1);
    error2=norm(w-w_esatta,2)/norm(w_esatta,2);
    errorinf=norm(w-w_esatta,'Inf')/norm(w_esatta,'Inf');
    V=V/norm(V);
    V_esatti=V_esatti/norm(V_esatti);
    e1(i)=error1;
    e2(i)=error2;
    einf(i)=errorinf;
end

%% Plot
cd ..
cd convergence_studies\Clamped_wing\ % Move to the right folder 

fig=figure;
set(gcf, 'Position',  [0, 0, 5000, 4000])
subplot(1,3,1)
    loglog(nel_tot,einf)
    grid on
    xlabel('Number of elements')
    ylabel('Norm INF error')
subplot(1,3,2)
    loglog(nel_tot,e2)
    grid on
    xlabel('Number of elements')
    ylabel('Norm 2 error')
subplot(1,3,3)
    loglog(nel_tot,t)
    grid on
    xlabel('Number of elements')
    ylabel('Time spent')
   
if 1
    saveas(fig,'Convergence_clamped_wing_eig','eps')
end
%% plot Modes
if 0
    options.plot_original = 1;
    options.plot_deformed = 1;
    options.plotColor = 'green';
    options.saveSTL = 0;
    options.point_section = 8;
    options.N = 30;
    m_Modes3d(model,options);
end