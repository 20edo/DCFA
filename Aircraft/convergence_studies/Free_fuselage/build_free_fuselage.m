function aircraft=build_free_fuselage(nel)
% Builds the fuselage given the number of elements

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


%% Define the nodes of the beam
aircraft=m_init();

Node1=en_free([0,0,0]);      % Where the node begins
Node2=en_free([6 0 0]);     % Where the section from conical becomes cylindrical
Node3=en_free([18 0 0]);    % Where the wings are clamped
Node4=en_free([40 0 0]);    % Where the conical section of the tail begins
Node5=en_free([48 0 2.5]);  % Where the rudder is clamped
Node6=en_free([52 0 2.5]);  % Where the aircraft ends

aircraft.en=[Node1, Node2, Node3, Node4, Node5, Node6];


%% Discretisation parameters
nel_nose=ceil(nel*Node2.x(1)/Node6.x(1));
nel_front_fuselage=round(nel*abs((Node3.x(1)-Node2.x(1))/Node6.x(1)));
nel_rear_fuselage=round(nel*abs((Node4.x(1)-Node3.x(1))/Node6.x(1)));
nel_front_tail=round(nel*abs((Node5.x(1)-Node4.x(1))/Node6.x(1)));
nel_rear_tail=nel-nel_nose-nel_front_fuselage-nel_rear_fuselage-nel_front_tail;

%% Parameters of the fuselage
R_fus= 3;       % Radius of the fuselage
t_fus=0.03;      % Thickness of the fuselage
R_mid_tail=1;   % Radius of the tail fuselage at the clamp with the rudder

%% Material properties (Aluminium)

E=70*1e9;   % Young modulus
G=27*1e9;   % Shear modulus
rho=2700;   % Density

%% Nose

L=norm(aircraft.en(1).x-aircraft.en(2).x,2);
R= @(x) R_fus*sqrt(1-((x-L)/L).^2);
t= @(x) R(x)/10;
nose=b_constant_p_tube(L,R,t,E,G,rho,nel_nose);
nose.o=aircraft.en(1).x;
nose.vx=[1 0 0]';
nose.vy=[0 1 0]';
nose.name='nose';

%% Front fuselage

L=norm(aircraft.en(2).x-aircraft.en(3).x,2);
R= @(x) R_fus+0.*x;
t= @(x) t_fus+0.*x;
front_fuselage=b_constant_p_tube(L,R,t,E,G,rho,nel_front_fuselage);
front_fuselage.o=aircraft.en(2).x;
front_fuselage.vx=[1 0 0]';
front_fuselage.vy=[0 1 0]';
front_fuselage.name='front_fuselage';

%% Rear fuselage

L=norm(aircraft.en(3).x-aircraft.en(4).x,2);
R= @(x) R_fus+0.*x;
t= @(x) t_fus+0.*x;
rear_fuselage=b_constant_p_tube(L,R,t,E,G,rho,nel_rear_fuselage);
rear_fuselage.o=aircraft.en(3).x;
rear_fuselage.vx=[1 0 0]';
rear_fuselage.vy=[0 1 0]';
rear_fuselage.name='rear_fuselage';

%% Front tail

L=norm(aircraft.en(4).x-aircraft.en(5).x,2);
R= @(x) R_fus-(R_fus-R_mid_tail)/L*x;
t= @(x) t_fus+0.*x;
front_tail=b_constant_p_tube(L,R,t,E,G,rho,nel_front_tail);
front_tail.o=aircraft.en(4).x;
front_tail.vx=aircraft.en(5).x-aircraft.en(4).x;
front_tail.vx=front_tail.vx/norm(front_tail.vx,2);
front_tail.vy=[0 1 0]';
front_tail.name='front_tail';

%% Rear tail

L=norm(aircraft.en(5).x-aircraft.en(6).x,2);
R= @(x)  R_mid_tail-(R_mid_tail-2*t_fus)/L*x+0.00001;
t= @(x) t_fus+0.*x;
rear_tail=b_constant_p_tube(L,R,t,E,G,rho,nel_rear_tail);
rear_tail.o=aircraft.en(5).x;
rear_tail.vx=aircraft.en(6).x-aircraft.en(5).x;
rear_tail.vx=rear_tail.vx/norm(rear_tail.vx,2);
rear_tail.vy=[0 1 0]';
rear_tail.name='rear_tail';

fuselage_beams=[nose, front_fuselage, rear_fuselage, front_tail, rear_tail];

%% Add beams to aircraft
for i=1:length(fuselage_beams)
    aircraft=m_add_beam(aircraft,fuselage_beams(i));
%     disp(fuselage_beams(i).name)
%     M=sum(sum(fuselage_beams(i).M(1:6:end,1:6:end)));
%     disp('M=')
%     disp(M)
    
end


%% Add payload
human=false;
vehicles=true;

% Cockpit and instruments

L=norm(aircraft.b(1).L);
l=1;
nel=(size(aircraft.b(1).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(1).o;
Payload.vx=aircraft.b(1).vx;
Payload.vy=aircraft.b(1).vy;
% disp('Cockpit')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(1).M=aircraft.b(1).M+Payload.M;

% Nose instuments

L=norm(aircraft.b(2).L);
l=1;
nel=(size(aircraft.b(2).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(2).o;
Payload.vx=aircraft.b(2).vx;
Payload.vy=aircraft.b(2).vy;
% disp('Nose ')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(2).M=aircraft.b(2).M+Payload.M;

if human
    
    % Material properties (Soldiers)

    E=0;            % Young modulus
    G=0;            % Shear modulus
    rho=1000/10;    % Density


    % Soldiers    
    L=norm(aircraft.b(3).L);
    l=2;
    nel=(size(aircraft.b(3).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(3).o;
    Payload.vx=aircraft.b(3).vx;
    Payload.vy=aircraft.b(3).vy;
%     disp('Soldiers_1')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(3).M=aircraft.b(3).M+Payload.M;
    
    
    % Soldiers    
    L=norm(aircraft.b(4).L);
    l=2;
    nel=(size(aircraft.b(4).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(4).o;
    Payload.vx=aircraft.b(4).vx;
    Payload.vy=aircraft.b(4).vy;
%     disp('Soldiers_2')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(4).M=aircraft.b(4).M+Payload.M;
    
    
elseif vehicles
    
    %% Material properties (Vehicles)

    E=0;            % Young modulus
    G=0;            % Shear modulus
    rho=8000/20;    % Density
    
    % Vehicles
    
    L=norm(aircraft.b(3).L);
    l=2.5;
    nel=(size(aircraft.b(3).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(3).o;
    Payload.vx=aircraft.b(3).vx;
    Payload.vy=aircraft.b(3).vy;
%     disp('Vehicles_1')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(3).M=aircraft.b(3).M+Payload.M;
    
    % Vehicles   
    L=norm(aircraft.b(4).L);
    l=2;
    nel=(size(aircraft.b(4).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(4).o;
    Payload.vx=aircraft.b(4).vx;
%     Payload.vy=aircraft.b(4).vy;
%     disp('Vehicles_2')
    disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(4).M=aircraft.b(4).M+Payload.M;
end

% APU and other systems

E=0;            % Young modulus
G=0;            % Shear modulus
rho=2700/20;    % Density

L=norm(aircraft.b(2).L);
l=1;
nel=(size(aircraft.b(2).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(2).o;
Payload.vx=aircraft.b(2).vx;
Payload.vy=aircraft.b(2).vy;
% disp('APU')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(2).M=aircraft.b(2).M+Payload.M;
