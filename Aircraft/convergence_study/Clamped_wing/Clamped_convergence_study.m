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
cd convergence_study\Clamped_wing\

%%

nel_tot=30:10:200;
t=zeros(size(nel_tot));
load w_esatta
number=20;      % Number of eigenvalues considered
w_esatta(number+1:end)=[];

for i =1:length(nel_tot)
    tic
    % Build model
    model=build_clamped_wing(aircraft,nel_tot(i));
    model=m_compute_matrices(model);
    % Solve
    [w, V, k] = ROM_solver(number, model.M, model.K);
    t(i)=toc;
    % Calculate errors
    error1=norm(w-w_esatta,1)/norm(w_esatta,1);
    error2=norm(w-w_esatta,2)/norm(w_esatta,2);
    errorinf=norm(w-w_esatta,'Inf')/norm(w_esatta,'Inf');
    e1(i)=error1;
    e2(i)=error2;
    einf(i)=errorinf;
    it(i)=k;
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
    plot(nel_tot,it)
    grid on
    title('Number of iterations')
    xlabel('Number of elements')
    ylabel('ROM_iterations')
   
saveas(fig,'Convergence_clamped_wing_eig','svg')