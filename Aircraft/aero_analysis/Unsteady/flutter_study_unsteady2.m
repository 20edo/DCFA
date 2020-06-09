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
n = 4;
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
v = [0:1:30];
q = 1/2*rho.*v.^2;


% First iteration
[X_old,e_old] = polyeig(K,Cs,M);
% I = imag(e_old)<=0.1;
% e_old = e_old.*I;
% e_pulito = [];
% X_pulito = [];
% for i = 1:2*n
%     if abs(e_old(i))>1e-3
%         e_pulito = [e_pulito; e_old(i)];
%         X_pulito = [X_pulito, X_old(:,i)];
%     end
% end
% X_old = X_pulito;
% e_old = e_pulito;
X_zero = X_old;
e_zero = e_old;

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

eig_ = zeros(length(v),size(X_old,2));
eig_(1,:) = e_old;


% Following iterations
for i=2:length(v)
    for k=1:length(e_old)
        tic
        %e_old(k)
        % Build matrices to solve the non-linear system
        kk = l*imag(e_old(k))/v(i);
        wing = m_add_unsteady_loads(wing,[1,0,0]',kk);
        Ham = wing.Ham;
        Ham = V'*Ham*V;
        [X,e] = polyeig(K-q(i)*Ham,Cs,M);
        %e
        d_e = abs(e-e_old(k));
        %d_e
        if 1 %sum(d_e < 0.1) == 1
            e = e(find(d_e == min(d_e)));
            X = X(:,find(d_e == min(d_e)));
        else
            guess = find(d_e < 0.1);
            d_v = vecnorm(X_old(:,guess)-X);
            index = guess(find(d_v == min(d_v)));
            e = e(index);
            X = X(:,index);
        end
        X_old(:,k) = X;
        e_old(k) = e;
        phrase = ['Eig number ',num2str(k),' out of ',num2str(length(e_old)),'; Velocity ',num2str(i),' out of ',num2str(length(v))];
        disp(phrase)
        toc
        
    end
    eig_(i,:) = e_old;
end

%% V-g plot
close all

g = 2*real(eig_)./(imag(eig_));
%% Plot frequency diagram and V-G diagram
figure
hold on
subplot(2,1,1)
plot(v,abs(imag(eig_)));
ylabel('imag(eig)')
grid on

subplot(2,1,2)
plot(v,g);
ylabel('g')
grid on
ylim([-0.05,0.05])


%% Plot the corresponding modeshapes
if 0
    figure
    phi = deg2rad(45);
    for i=4:5
        wing.b(i).ssh = true;
    end
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    
    num = sum(imag(e_zero)>=-0.1);
    i = 1;
    for k = 1:2*n
        if imag(e_zero(k))>=-0.1
            subplot(2,(num+mod(num,2))/2,i)
            m_plot_eigenshape2(wing,options,real(exp(1i*phi)*V*(X_zero(:,k))*30));
            title(num2str(abs(imag(e_zero(k)))))
            i = i+1;
        end
    end
end