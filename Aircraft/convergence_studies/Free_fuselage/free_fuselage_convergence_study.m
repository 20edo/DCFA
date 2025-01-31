% This script is made to study the convergence of eigenvalues of the
% free-free fuselage

%%
clear all, close all, clc

%% Generate the aircraft model
cd ..
cd ..
cd generate_model\
generate_model

%% Load the exact solution

cd convergence_studies\Free_fuselage\ % Move to the right folder 

load w_esatta.mat
load V_esatti.mat

cd ..
cd ..
cd generate_model

%% Perform analysis
save=1;
% nel_tot=500;          % Exact (reference) eigenvalues
nel_tot=50:10:300;
t=zeros(size(nel_tot));
number=100;      % Number of eigenvalues considered
w_esatta(number+1:end)=[];
V_esatti(:,number+1:end)=[];
model=m_init();

for i =1:length(nel_tot)
    tic
    % Build model
    model=build_free_fuselage(nel_tot(i));
    model=m_compute_matrices(model);
    % Shift
    alpha=1;
    % Solve
    [V,D,flag] = eigs(model.K+alpha*model.M,model.M,number+6,'smallestabs');
    w=real(diag(D-alpha).^0.5);
    t(i)=toc;
    % Suppress rigid body egienvalues
    w(1:6)=[];
    % Calculate errors
    error1=norm(w-w_esatta,1)/norm(w_esatta,1);
    error2=norm(w-w_esatta,2)/norm(w_esatta,2);
    errorinf=norm(w-w_esatta,'Inf')/norm(w_esatta,'Inf');
    e1(i)=error1;
    e2(i)=error2;
    einf(i)=errorinf;
end

cd ..
cd convergence_studies\Free_fuselage\ % Move to the right folder 

%% Plot

fig=figure;
set(gcf, 'Position', [0, 0, 600, 400])
loglog(nel_tot,einf)
grid on
xlabel('Number of elements','interpreter','latex')
ylabel('Norm INF error','interpreter','latex')
if save
    saveas(fig,'Norminf_convergence_free_fuselage','epsc')
end
fig=figure;
set(gcf, 'Position', [0, 0, 600, 400])
loglog(nel_tot,e2)
grid on
xlabel('Number of elements','interpreter','latex')
ylabel('Norm 2 error','interpreter','latex')
if save
    saveas(fig,'Norm2_convergence_free_fuselage','epsc')
end
fig=figure;
set(gcf, 'Position', [0, 0, 600, 400])
loglog(nel_tot,t)
grid on
xlabel('Number of elements','interpreter','latex')
ylabel('Time spent \quad [s]','interpreter','latex')
if save
    saveas(fig,'Time_convergence_free_fuselage','epsc')
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