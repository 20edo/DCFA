% Consistent problems


clear all, close all, clc
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
for i=1:length(wing_list)
    wing=m_add_beam(wing,wing_list(i));
end

%% Add the aero loads
wing = m_add_aero_loads(wing,[1,0,0]');

%% If you want to study different position of the ailerons
wing.b(1).fb = zeros(size(wing.b(1).fb));
wing.b(1).Lq = zeros(size(wing.b(1).Lq));
wing.b(1).Lb = zeros(size(wing.b(1).Lb));
wing.b(1).lb = zeros(size(wing.b(1).lb));
wing.b(2).fb = zeros(size(wing.b(2).fb));
wing.b(2).Lq = zeros(size(wing.b(2).Lq));
wing.b(2).Lb = zeros(size(wing.b(2).Lb));
wing.b(2).lb = zeros(size(wing.b(2).lb));

wing = m_compute_matrices(wing);

%% Compute J_x of the aircraft
Jx = 0;         % Inizialize Jx
% External nodes
for i=1:length(aircraft.en)
   Jx=Jx+aircraft.en(i).M(4,4)+aircraft.en(i).M(1,1)*norm(aircraft.en(1).x([2,3]))^2; 
end

% Beams
for i=1:length(aircraft.b)
    beam=aircraft.b(i);
    for j=0:length(beam.in)-1
        x=beam.o+beam.vx*beam.in(j+1).x;
        Jx=Jx+beam.M(4+6*j,4+6*j)+beam.M(1+6*j,1+6*j)*norm(x([2,3]))^2;
    end
end

lp = wing.lp;
lb = wing.lb;
lq = wing.lq;
Sq = wing.Sq;
fp = wing.fp;
fb = wing.fb;
K = wing.K;
Ka = wing.Ka;

%% Number of problem to be solved
problem=3;
%% #1
% initial roll acceleration p_dot for prescribed aileron deflection beta
%% #2
% aileron deflection b, for prescribed initial roll acceleration p_dot
%% #3
% asymptotic roll rate for prescribed beta
%% #4 
% aileron deflection beta for prescribed asymptotic roll rate 
%% #5 
% initial roll acc p_dot for prescribed roll rate p
%% 6 
% roll rate p required for prescribed roll acceleration 

q = 24500; % dynamic pressure
switch problem
    case 1
        beta = deg2rad(2);
        A = [K-q*Ka, Sq
            -q*lq, Jx];
        b = q*[fb; lb]*beta;
        sol = A\b;
        p_dot = sol(end);
        [V,D] = eig(full([K, Sq; zeros(1,size(K,1)), Jx]),full([Ka, zeros(size(Ka,1),1); lq, 0]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the divergence dynamic pressure for the manouver')
        end
        
    case 2
        p_dot = 1;
        A = [K-q*Ka, -q*fb; -q*lq, -q*lb];
        b = -[Sq; Jx]*p_dot;
        sol = A\b;
        beta = sol(end);
        [V,D] = eig(full([K, zeros(size(K,1),1); zeros(1,size(K,1)), 0]),full([Ka, fb; lq, lb]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the control reversal dynamic pressure for the manouver')
        end
        
    case 3
        beta = deg2rad(20);
        A = [K-q*Ka, -q*fp; -q*lq, -q*lp];
        b = q*[fb; lb]*beta;
        sol = A\b;
        p_fract_vinf = sol(end);
        [V,D] = eig(full([K, zeros(size(K,1),1); zeros(1,size(K,1)), 0]),full([Ka, fp; lq, lp]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the divergence dynamic pressure for the manouver')
        end
    case 4 
        p_fract_vinf = 0.05; 
        A = [K-q*Ka, -q*fb; -q*lq, -q*lb]; 
        b = q*[fp; lp]*p_fract_vinf; 
        sol = A\b; 
        beta = sol(end); 
        [V,D] = eig(full([K, zeros(size(K,1),1); zeros(1,size(K,1)), 0]),full([Ka, fb; lq, lb]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the control reversal dynamic pressure for the manouver')
        end
    case 5 
        p_fract_vinf = 0.05; 
        A = [K-q*Ka, Sq; -q*lq, Jx]; 
        b = q*[fp; lp]*p_fract_vinf; 
        sol = A\b; 
        p_dot = sol(end); 
        [V,D] = eig(full([K, Sq; zeros(1,size(K,1)), Jx]),full([Ka, zeros(size(Ka,1),1); lq, 0]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the divergence dynamic pressure for the manouver')
        end
    case 6 
        p_dot = 0.04; 
        A = [K-q*Ka, -q*fp; -q*lq, -q*lp]; 
        b = -[Sq; Jx]*p_dot; 
        sol = A\b; 
        p_fract_vinf = sol(end); 
        [V,D] = eig(full([K, zeros(size(K,1),1); zeros(1,size(K,1)), 0]),full([Ka, fp; lq, lp]));
        q_div = diag(D);
        [q_div,I] = sort(real(q_div));
        q_div(q_div<1)=[];
        q_div = q_div(1);
        if q > q_div
            warning('You are over the divergence dynamic pressure for the manouver')
        end     
end

%% Control aeroelastic

q_ref=linspace(1,15000*1e3,5*1e2);
elastic_control=zeros(size(q_ref));

for i=1:length(q_ref)
    % Elastic solution, problem 1
    q=q_ref(i);
    beta = deg2rad(2);
    A = [K-q*Ka, Sq
        -q*lq, Jx];
    b = q*[fb; lb]*beta;
    sol = A\b;
    p_dot(i) = sol(end);
    
    % Rigid solution problem 1
    
    A_r = [K, Sq
           -0*lq, Jx];
    b_r = q*[fb; lb]*beta;
    sol_r = A_r\b_r;
    p_dot_r(i) = sol_r(end);
    elastic_control(i)=p_dot(i)/p_dot_r(i);
    
    
end

% %% Roll damping aeroelastic correction
% 
% elastic_damping=zeros(size(q_ref));
% 
% for i=1:length(q_ref)
%     q=q_ref(i);
%     
%     % Elastic solution, problem 5
%     p_fract_vinf = deg2rad(2)/200;
%     A = [K-q*Ka, Sq; -q*lq, Jx];
%     b = q*[fp; lp]*p_fract_vinf;
%     sol = A\b;
%     p_dot = sol(end);
%     
%     % Rigid solution, problem 5
%     
%     A_r = [K, Sq
%         -0*lq, Jx];
%     b_r= q*[fp; lp]*p_fract_vinf;
%     sol_r = A_r\b_r;
%     p_dot_r = sol_r(end);
%     elastic_damping(i)=p_dot/p_dot_r;
% end


%% Plot

figure 
hold on 
plot(q_ref,p_dot,'b')
plot(q_ref,p_dot_r,'r')
hold off 
 
% Deformative shape plot

if 0
    % switch on the aero properties for the plot 
    for i=4:5
        wing.b(i).ssh = true; 
    end
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 1;        % we have only one eig
    for g = 1:length(sol)-1
        if mod(g,6)==3 || mod(g,6)==5
            sol(g) = -sol(g); 
        end
    end
    m_plot_eigenshape(wing,options,sol(1:end-1)*40)
end
    
% Control aeroelastic correction

if 1
    figure 
    plot(q_ref,elastic_control)
    xlabel('Dynamic pressure')
    ylabel('$\frac{\dot{p}}{\dot{p_r}}$','Interpreter','latex')
    title('Roll control aeroelastic correction (steady aerodynamics)')
    grid on
    ylim([-5 30])
end

% % Damping aeroelastic correction
% if 1
%     figure 
%     plot(q_ref,elastic_damping)
%     xlabel('Dynamic pressure')
%     ylabel('$\frac{\dot{p}}{\dot{p_r}}$','Interpreter','latex')
%     title('Roll aeroelastic damping correction (steady aerodynamics)')
%     grid on
% end








