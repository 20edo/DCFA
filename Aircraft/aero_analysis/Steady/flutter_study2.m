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
[nn,II] = sort(imag(e_old));
e_old = e_old(II);
X_old = X_old(:,II);
e_old = e_old(length(e_old)/2+1:end);
X_old = X_old(:,size(X_old,2)/2+1:end);
X_zero = X_old;
e_zero = e_old;

% Initialize non linear system variables
A=zeros(size(M,1)+1);
b=zeros(size(M,1)+1,1);

eig_ = zeros(length(v),length(e_old));
eig_(1,:) = e_old;
X_save = zeros(length(v),size(X_old,1),size(X_old,2));
X_save(1,:,:) = X_old;
% Following iterations
for i=2:length(v)
    X = zeros(size(X_old));
    e = zeros(size(e_old));
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
%% Plot frequency diagram and V-G diagram all togheter
if 0
    figure(1)
    hold on
    subplot(3,1,1)
    plot(v,abs(imag(eig_))/(2*pi),'LineWidth',1.5);
    ylabel('Frequency \quad [Hz]','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',18,'interpreter','latex');
    grid on
    
    subplot(3,1,2)
    plot(v,g,'LineWidth',1.5);
    ylabel('g','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    grid on
    ylim([-0.3,0.15])
    
    subplot(3,1,3)
    plot(v,g,'LineWidth',1.5);
    ylabel('g','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    grid on
    ylim([-0.06,0.04])
    xlim([0,400])
    set(gcf, 'Position',  [0, 0, 700, 800])
%     saveas(figure(1),'flutter_1','epsc')
end
%% V-G diagram separated
if 0
    figure(1)
    hold on
    plot(v,abs(imag(eig_))/(2*pi),'LineWidth',1.5);
    ylabel('Frequency \quad [Hz]','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    grid on
    set(gcf, 'Position',  [0, 0, 700, 250])
    saveas(figure(1),'flutter_5','epsc')
    
    figure(2)
    hold on
    plot(v,g,'LineWidth',1.5);
    ylabel('g','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    grid on
    ylim([-0.3,0.15])
    set(gcf, 'Position',  [0, 0, 700, 250])
    saveas(figure(2),'flutter_6','epsc')
    
    figure(3)
    plot(v,g,'LineWidth',1.5);
    ylabel('g','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    grid on
    ylim([-0.06,0.04])
    xlim([0,400])
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    set(gcf, 'Position',  [0, 0, 700, 250])
    saveas(figure(3),'flutter_7','epsc')
end
%% Plot the partecipation of each mode 
if 0
    figure(4)
    plot(v,abs(X_save(:,:,3)),'LineWidth',1.5)
    ylabel('Mode contribution','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('$3^{rd}$ mode','fontsize',14,'interpreter','latex');
    grid on
    hold on 
    p = plot(286*ones(1e3,1),linspace(0,1.2,1e3),'Color','k');
    legend(p(1),'Flutter velocity')
%     xlim([2,1000])
    set(gcf, 'Position',  [0, 0, 500, 400])
%     saveas(figure(4),'flutter_8','epsc')
    
 
    figure(5)
    plot(v,abs(X_save(:,:,4)),'LineWidth',1.5)
    ylabel('Mode contribution','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('$4^{rd}$ mode','fontsize',14,'interpreter','latex');
    grid on
    hold on 
    p = plot(277*ones(1e3,1),linspace(0,1.2,1e3),'Color','k');
    legend(p(1),'Flutter velocity')
%     xlim([2,1000])
    set(gcf, 'Position',  [0, 0, 500, 400])
%     saveas(figure(5),'flutter_9','epsc')
end
%% Plot the value of the reduced frequency
if 0
    chord=7.72;
    l = chord/2;
    red_freq = zeros(size(g));
    for i=1:size(eig_,2)
        red_freq(:, i) = l*abs(imag(eig_(:,i)))./v';
    end
    figure(2)
    plot(v,red_freq,'LineWidth',1.5)
    ylabel('k \quad [-]','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    ylim([0,0.5])
    grid on
    set(gcf, 'Position',  [40, 40, 500, 500])
    saveas(figure(2),'flutter_4','epsc')
end
%% Plot the corresponding modeshapes
if 0
    close all
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
            freq = [num2str(abs(imag(e_zero(k)))/(2*pi)),' [Hz]'];
            m_plot_eigenshape2(wing,options,real(exp(1i*phi)*V*(X_zero(:,k))*30));
            view(-65,35)
            title(freq,'fontsize',18,'interpreter','latex')
            i = i+1;
        end
    end
    set(gcf, 'Position',  [0, 0, 2000, 2000])
end

%% Plot the flutter eigenshape
for i=2:length(v)
    if  ~isempty(find(g(i,:)>0))
        break
    end
end
k = find(g(i,:)>0);
k = k(1);
flutter_mode = X_save(i,:,k)';
disp('flutter velocity')
disp(v(i));
disp('flutter mach')
disp(v(i)/a)

if 0
    angle = 40;
    if 0
        angle = 0:5:360;
    end
    for g = 1:length(angle)
        disp(angle(g))
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
        set(gcf, 'Position',  [0, 0, 2000, 2000])
        FR(g) = getframe(figure(1));
    end
    %%
    if 0
        % create the video writer with 1 fps
        writerObj = VideoWriter('myVideo3.avi','Uncompressed AVI');
        writerObj.FrameRate = 100;
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
    end
end

