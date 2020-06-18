% This script builds the typical section matrices for the unsteady
% aerodynamics.

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
clear all,close all, clc

%Add to path the necessary folders
cd ..
cd ..
cd ..
init
cd Aircraft\aero_analysis\Unsteady

%% Generate model of a section clamped for all gdl except rotation around x and torsion

% Material proerties
E=70*1e9;
G=27*1e9;
rho=2700;
% Edit to achieve insane stiffness
E=E*1e8;
G=G*1e8;
rho=rho;

% Geometrical properties of the section
c=@(x) 1+0.*x;
h=@(x) 0.12+0.*x;
t=@(x) 0.001+0.*x; % Has no sense but I want it very stiff-rigid)
L=1;
nel=10;

% Build beam
beam=b_ssh_profile(L,c,h,t,E,G,rho,nel);
beam.vx=[0 -1 0]';
beam.vy=[1 0 0]';

% Clamp the root (except some gdl)
node=en_ground([0,0,0]');
node.c(3)=false;
node.c(4)=false;
node.K=eye(6)*1e3;


model=m_init();
m.en=node;
model=m_add_beam(model,beam);

wing=m_compute_matrices(model);

%% Same as study fsolve
chord=1;
l = chord/2;

%% Reduction of the model using n eigenvectors
n = 6;
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
v = [100:10:200];
q = 1/2*rho.*v.^2;


scaling=1;  %redifined later

% First iteration
[X_old,e_old] = polyeig(K/scaling,Cs/scaling,M/scaling);
e_old=e_old*scaling;
I = imag(e_old)>-00000000.1;
e_old = e_old.*I;
e_pulito = [];
X_pulito = [];
for i = 1:2*n
    if abs(e_old(i))>1e-3
        e_pulito = [e_pulito; e_old(i)];
        X_pulito = [X_pulito, X_old(:,i)];
    end
end
X_old = X_pulito;
e_old = e_pulito;
[nn,II] = sort(imag(e_old));
e_old = e_old(II); 
X_old = X_old(:,II);
X_zero = X_old;
e_zero = e_old;
% scaling=sum(abs(e_zero));

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
exitflag=zeros(length(v),length(e_old));
t=zeros(length(v),length(e_old));

X_save = zeros(length(v),size(X_old,1),size(X_old,2));
% Following iterations
for i=2:length(v)
    for k=1:length(e_old)
        tic
        options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-8,'algorithm','levenberg-marquardt',...
            'ScaleProblem','jacobian','parallel',true);
        [x,~,exitflag(i,k)] = fsolve(@(Unknown) funz(wing,v(i),V,q(i),scaling,Unknown(2:end),Unknown(1)),[e_old(k);X_old(:,k)],options);
        e=x(1);
        X=x(2:end);
%         if 1 %sum(d_e < 0.1) == 1
%             e = e(find(d_e == min(d_e)));
%             X = X(:,find(d_e == min(d_e)));
%         else
%             guess = find(d_e < 0.1);
%             d_v = vecnorm(X_old(:,guess)-X);
%             index = guess(find(d_v == min(d_v)));
%             e = e(index);
%             X = X(:,index);
%         end
        X_old(:,k) = X;
        e_old(k) = e*scaling;
        phrase = ['Eig number ',num2str(k),' out of ',num2str(length(e_old)),'; Velocity ',num2str(i),' out of ',num2str(length(v))];
        disp(phrase)
        t(i,k)=toc;
    end
    eig_(i,:) = e_old;
    X_save(i,:,k) = X; 
%     scaling=sum(abs(eig_(i,:)));
end

%% V-g plot
close all

g = 2*real(eig_)./abs(imag(eig_));

% figure
% hold on
% for j = 2:2:size(eig_,2)
%     
% subplot(2,1,1)
% hold on 
% plot(v,abs(imag(eig_(:,j))))
% hold off 
% ylabel('imag(eig)')
% grid on
% 
% subplot(2,1,2)
% hold on 
% plot(v,g(:,j))
% hold off 
% ylabel('g')
% grid on
% end

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
if 1
    figure
    phi = deg2rad(45);
    clear options
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

load handel
sound(y,Fs)