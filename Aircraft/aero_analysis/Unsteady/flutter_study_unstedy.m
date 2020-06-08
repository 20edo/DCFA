% This is a script to study the flutter problem at 10.000 m for the clamped
% wing

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%
%
%
clear all , close all, clc
cd ..
cd ..
%% Generate the wing model
cd generate_model
generate_model
% Move to the analysis folder
cd aero_analysis\Unsteady

% switch off the aerodynamic properties of the engine support
for i=16:19
    aircraft.b(i).ssh = false;
end

%% Build the swept wing model
wing=m_init();
wing.en=[en_ground(aircraft.en(7).x) ...
    aircraft.en(17) aircraft.en(18)];
wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9) aircraft.b(16) aircraft.b(17)];
% wing.en=[en_ground(aircraft.en(7).x)];
% wing_list=[aircraft.b(7) aircraft.b(8) aircraft.b(9)];
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end

chord=7.72;
l = chord/2;

wing = m_compute_matrices(wing);
%% Reduction of the model using n eigenvectors
n = 18;
[V,D] = eigs(wing.K,wing.M,n,'smallestabs');
V_red = V;

alpha = 0;
gamma = 0;
Cs = alpha*wing.M + gamma*wing.K;
% Cs = 1e-3*sum(sum(diag(wing.K)))/size(wing.K,1)*ones(size(wing.K));  

M = V'*wing.M*V;
K = V'*wing.K*V;
Cs = V'*Cs*V; 

%% Altitude fixed to 10.000 m
[T,a,P,rho] = atmosisa(10000);

% the problem is in the form
% M*q_dotdot - q/Vinf*C*q_dot + (K - q*Ka)*q = 0
% the solution of the problem is given by polyeig(K,C,M)

%% Tracking of eigenvalues trough eigenvectors
v = [0:15:400];
q = 1/2*rho.*v.^2;


% First iteration
[X_old,e_old] = polyeig(K,Cs,M);

% % Find the derivatives of Ham
% k1 = 0.5e-12;
% wing = m_add_unsteady_loads(wing,[1,0,0]',k1);
% Ham = wing.Ham;
% wing = m_add_unsteady_loads(wing,[1,0,0]',0);
% Ham_zero = wing.Ham;
% Ham_dk = 1i*imag(Ham)/k1;
% Ham_dk2 = 2*(real(Ham)-Ham_zero)/k1^2;
% 
% % Reduce matrices
% Ham_zero = V'*Ham_zero*V;
% Ham_dk = V'*Ham_dk*V;
% Ham_dk2 = V'*Ham_dk2*V;



% Initialize non linear system variables
A=zeros(size(M,1)+1);
b=zeros(size(M,1)+1,1);

eig_ = zeros(length(v),2*n); 
eig_(1,:) = e_old;


% Following iterations
for i=2:length(v)
    X = zeros(size(M,1),2*size(M,1));
    e = zeros(2*size(M,1),1);
    for k=1:size(X_old,2)
        tic
        kk = l*imag(e_old(k))/v(i); 
        wing = m_add_unsteady_loads(wing,[1,0,0]',kk); 
        Ham = wing.Ham;
        Ham_dk = wing.Ham_dk;
        Ham = V'*Ham*V; 
        Ham_dk = V'*Ham_dk*V;
        A(1:size(M,1),1:size(M,1))=e_old(k)^2*M+e_old(k)*Cs+K-q(i)*(Ham);
        A(1:size(M,1),end)=(2*e_old(k)*M+Cs-q(i)*(-1i*Ham_dk*l/v(i)))*X_old(:,k);
        A(end,1:size(M,1))=2*X_old(:,k)';
        A(end,end)=0;
        b(1:size(M,1),1)=-A(1:size(M,1),1:size(M,1))*X_old(:,k);
        b(end)=1-X_old(:,k)'*X_old(:,k);
%         funz=@(z) A*z-b;
%         z0=[X_old(:,k);e_old(k)];
%         [z,~,exitflag]=fsolve(funz,z0);
        z=A\b;
        X(:,k)=X_old(:,k)+z(1:end-1);
        e(k)=e_old(k)+z(end);
        phrase = ['Eig number ',num2str(k),' out of ',num2str(2*n),'; Velocity ',num2str(i),' out of ',num2str(length(v))];
        disp(phrase)
        toc
        
    end
    X_old = X;
    e_old = e;
    eig_(i,:) = e_old;
end

%% V-g plot
close all 
g = 2*real(eig_)./abs(imag(eig_));
figure
hold on
ylim([-0.5,0.5])

plot(v,g);
    
ylim([-4,0.5])
ylabel('g')
grid on 
title('V-G diagram')

% figure 
% hold on 
% for k = 1:n
%     plot(v,imag(eig_(:,k)));
% end
% ylabel('w')
% grid on 
% title('frequency-V diagram')



% figure
% hold on
% for k = 1:n
%     plot(v,real(eig_(:,k)));
% end
% ylim([-50,50])
% ylabel('real(eig)')


for i=4:5
        wing.b(i).ssh = true; 
    end
    
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    m_plot_eigenshape(wing,options,real(V*X(:,31)*30));
%     m_plot_eigenshape(wing,options,real(V*X(:,28))*30);



