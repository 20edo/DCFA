% This script is made to study the convergence of eigenvalues of the
% free-free fuselage

%%

nel_tot=300;
t=zeros(size(nel_tot));
% load w_esatta
w_esatta=zeros(50,1);
number=50;      % Number of eigenvalues considered
w_esatta(number+1:end)=[];

for i =1:length(nel_tot)
    tic
    % Build model
    model=build_free_fuselage(nel_tot(i));
    model=m_compute_matrices(model);
    % Shift
    alpha=1;
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
    it(i)=k;
end

%% Plot
fig=figure;

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
    ylabel('ROM iterations')
   
saveas(fig,'Convergence_free_fuselage_eig','svg')