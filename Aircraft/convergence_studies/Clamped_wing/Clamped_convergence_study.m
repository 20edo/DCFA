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

%% Generate the aircraft model
cd ..
cd ..
cd generate_model\
generate_model
cd convergence_studies\Clamped_wing\

%%

nel_tot=400;
t=zeros(size(nel_tot));
load w_esatta.mat
load V_esatti.mat
number=100;      % Number of eigenvalues considered
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
    e1(i)=error1;
    e2(i)=error2;
    einf(i)=errorinf;
    e_vect(i)=norm(V(24:30,:)-V_esatti(24:30,:));
end

%% Plot
fig=figure

subplot(2,2,1)
    loglog(nel_tot,einf)
    grid on
    xlabel('Number of elements')
    ylabel('Error')
    title('Norm INF error')
subplot(2,2,2)
    loglog(nel_tot,e2)
    grid on
    xlabel('Number of elements')
    ylabel('Error')
    title('Norm 2 error')
subplot(2,2,3)
    loglog(nel_tot,t)
    grid on
    xlabel('Number of elements')
    ylabel('Time')
    title('Time spent')
subplot(2,2,4)
    plot(nel_tot,e_vect)
    grid on
    title('Norm 2 tip displacements error')
    xlabel('Number of elements')
    ylabel('Norm 2 error')
   
saveas(fig,'Convergence_clamped_wing_eig','svg')