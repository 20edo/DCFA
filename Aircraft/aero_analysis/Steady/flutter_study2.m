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
cd aero_analysis\Steady

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
%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');

%% Reduction of the model using n eigenvectors
n = 6;
[V,D] = eigs(wing.K,wing.M,n,'smallestabs');
V_red = V;

alpha = 0;
gamma = 0;
Cs = alpha*wing.M + gamma*wing.K;
% Cs = 1e-3*sum(sum(diag(wing.K)))/size(wing.K,1)*eye(size(wing.M));

M = V'*wing.M*V;
K = V'*wing.K*V;
Ka = V'*wing.Ka*V;
Ca = V'*wing.Ca*V;
Cs = V'*Cs*V;

%% Altitude fixed to 10.000 m
[T,a,P,rho] = atmosisa(10000);

% the problem is in the form
% M*q_dotdot - q/Vinf*C*q_dot + (K - q*Ka)*q = 0
% the solution of the problem is given by polyeig(K,C,M)

%% Tracking of eigenvalues trough eigenvectors
v = [0:1:1000];
q = 1/2*rho.*v.^2;

% Cs = 1e-3*sum(sum(diag(K)))/size(K,1)*eye(size(M));
% % alpha = 0.1;
% % gamma = 0.1;
% % Cs = alpha*M + gamma*K;

% First iteration

[X_old,e_old] = polyeig(K,Cs,M);
X_zero = X_old;
e_zero = e_old;
% Initialize non linear system variables
A=zeros(size(M,1)+1);
b=zeros(size(M,1)+1,1);

eig_ = zeros(length(v),2*n);
eig_(1,:) = e_old;
X_save = zeros(length(v),size(X_old,1),size(X_old,2));

% Following iterations
for i=2:length(v)
    X = zeros(size(M,1),2*size(M,1));
    e = zeros(2*size(M,1),1);
    for k=1:size(X_old,2)
        A(1:size(M,1),1:size(M,1))=e_old(k)^2*M+e_old(k)*(Cs-q(i)/v(i)*Ca)+K-q(i)*Ka;
        A(1:size(M,1),end)=(2*e_old(k)*M+Cs-q(i)/v(i)*Ca)*X_old(:,k);
        A(end,1:size(M,1))=2*X_old(:,k)';
        A(end,end)=0;
        b(1:size(M,1),1)=-A(1:size(M,1),1:size(M,1))*X_old(:,k);
        b(end)=1-X_old(:,k)'*X_old(:,k);
        z=A\b;
        X(:,k)=X_old(:,k)+z(1:end-1);
        e(k)=e_old(k)+z(end);
    end
    
    X_old = X;
    e_old = e;
    eig_(i,:) = e_old;
    X_save(i,:,:) = X; 
end

%% V-g plot
close all

g = 2*real(eig_)./abs(imag(eig_));
%% Plot frequency diagram and V-G diagram
figure(1)
hold on
subplot(3,1,1)
plot(v,abs(imag(eig_))/(2*pi));
ylabel('Frequency \quad [Hz]','fontsize',14,'interpreter','latex')
xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
grid on

subplot(3,1,2)
plot(v,g);
ylabel('g','fontsize',14,'interpreter','latex')
xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
grid on
ylim([-0.3,0.15])

subplot(3,1,3)
plot(v,g);
ylabel('g','fontsize',14,'interpreter','latex')
xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
grid on
ylim([-0.06,0.04])
set(gcf, 'Position',  [0, 0, 700, 800])
saveas(figure(1),'flutter_1','epsc')

%% Plot the corresponding modeshapes
if 0
    figure
    phi = deg2rad(50);
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
    
    [~, II] = sort(imag(e_zero));
    e_zero = e_zero(II);
    X_zero = X_zero(:,II);
    
    for k = 1:2*n
        if imag(e_zero(k))>=-0.1
            subplot(2,(num+mod(num,2))/2,i)
            freq = [num2str(abs(imag(e_zero(k)))/(2*pi)),' Hz']; 
            m_plot_eigenshape2(wing,options,real(exp(1i*phi)*V*(X_zero(:,k))*30));
            title(freq)
            i = i+1;
        end
    end
    
end

%% Plot the flutter eigenshape
for i=1:length(v)
   if  ~isempty(find(g(i,:)>0))
       break
   end
end
k = find(g(i,:)>0); 
k = k(1); 
flutter_mode = X_save(i,:,k)';

angle = 0:5:360;
for g = 1:length(angle)
close all 
    phi = deg2rad(angle(g));
    for i=4:5
        wing.b(i).ssh = true;
    end
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
m_plot_eigenshape2(wing,options,real(exp(1i*phi)*V*flutter_mode*15));
figure(1)
zlim([-10,10]);
xlim([10,40]); 
ylim([-5,35]);
view(-75,35)
FR(g) = getframe(figure(1));
end
%% 
  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
F = [FR,FR,FR,FR,FR,FR]; 
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
