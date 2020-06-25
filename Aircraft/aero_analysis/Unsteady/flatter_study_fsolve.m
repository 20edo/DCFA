% This is a script to study the unsteady flutter problem at 10.000 m for the clamped
% wing using the built in function of MATLAB

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
wing = m_add_aero_loads(wing,[1,0,0]');

chord=7.72;
l = chord/2;

wing = m_compute_matrices(wing);
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
v = [0:15:800];
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
X_save(1,:,:) = X_old;
% Following iterations
for i=2:length(v)
    for k=1:length(e_old)
        tic
        options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-8,'algorithm','levenberg-marquardt',...
            'ScaleProblem','jacobian','UseParallel',true);
        [x,~,exitflag(i,k)] = fsolve(@(Unknown) funz(wing,v(i),V,q(i),scaling,Unknown(2:end),Unknown(1)),[e_old(k);X_old(:,k)],options);
        e=x(1);
        X=x(2:end)/norm(x(2:end));
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
    X_save(i,:,:) = X_old;
    %     scaling=sum(abs(eig_(i,:)));
end

%% V-g plot all togheter
if 0
    close all
    
    g = 2*real(eig_)./abs(imag(eig_));
       
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
end
%% V-G diagram separated
if 1
    g = 2*real(eig_)./abs(imag(eig_));
    figure(1)
    hold on
    plot(v,abs(imag(eig_))/(2*pi),'LineWidth',1.5);
    ylabel('Frequency \quad [Hz]','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    grid on
    set(gcf, 'Position',  [0, 0, 700, 250])
%     saveas(figure(1),'un_flutter_12','epsc')
    
%     figure(2)
%     hold on
%     plot(v,abs(imag(eig_))/(2*pi),'LineWidth',1.5);
%     ylabel('Frequency \quad [Hz]','fontsize',14,'interpreter','latex')
%     xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
%     title('h = $10000$ m','fontsize',14,'interpreter','latex');
%     grid on
%     ylim([0,0.75])
%     set(gcf, 'Position',  [0, 0, 700, 250])
%     saveas(figure(2),'un_flutter_8','epsc')
    
    figure(3)
    hold on
    plot(v,g,'LineWidth',1.5);
    ylabel('g','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    grid on
    ylim([-0.25,0.15])
    set(gcf, 'Position',  [0, 0, 700, 250])
%     saveas(figure(3),'un_flutter_13','epsc')
    
%     figure(4)
%     plot(v,g,'LineWidth',1.5);
%     ylabel('g','fontsize',14,'interpreter','latex')
%     xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
%     grid on
%     ylim([-0.04,0.02])
%     title('h = $10000$ m','fontsize',14,'interpreter','latex');
%     set(gcf, 'Position',  [0, 0, 700, 250])
%     saveas(figure(4),'un_flutter_10','epsc')
end



%% Plot the partecipation of each mode 
if 1
    figure(5)
    plot(v,abs(X_save(:,:,3)),'LineWidth',1.5)
    ylabel('Mode contribution','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('$3^{rd}$ mode','fontsize',14,'interpreter','latex');
    grid on
    hold on
    p = plot(300*ones(1e3,1),linspace(0,1.2,1e3),'Color','k');
    legend(p(1),'Flutter velocity','Location','northeast')
%     xlim([2,1000])
    set(gcf, 'Position',  [0, 0, 500, 400])
    saveas(figure(5),'un_flutter_14','epsc')
end

%% Plot the reduced frequences 
if 0
    chord=7.72;
    l = chord/2;
    red_freq = zeros(size(g));
    for i=1:size(eig_,2)
        red_freq(:, i) = l*abs(imag(eig_(:,i)))./v';
    end
    figure(6)
    plot(v,red_freq,'LineWidth',1.5)
    ylabel('k \quad [-]','fontsize',14,'interpreter','latex')
    xlabel('VTAS \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
    title('h = $10000$ m','fontsize',14,'interpreter','latex');
    ylim([0,0.4])
    grid on
    set(gcf, 'Position',  [40, 40, 500, 500])
    saveas(figure(6),'un_flutter_6','epsc')
end

%% Plot the corresponding modeshapes
if 0
    figure
    phi = deg2rad(45);
    for i=4:5
        wing.b(i).ssh = true;
    end
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

%% Plot the flutter eigenshape
clearvars -except X_save g eig_ v a V wing
if 1
    for i=1:length(v)
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
end
%%
if 1
    angle = 145;
    if 1
        angle = 0:5:360;
    end
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
        set(gcf, 'Position',  [0, 0, 2000, 2000])
        FR(g) = getframe(figure(1));
    end
    %%
    if 1
        % create the video writer with 1 fps
        writerObj = VideoWriter('myVideo3.avi');
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
   


% load handel
% sound(y,Fs)