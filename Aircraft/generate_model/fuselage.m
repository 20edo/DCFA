% This script builds the beam that represent the fuselage of the aircraft
% The origin is set at the nose of the aircraft, the x axis is aligned with
% the fuselage and points to the front

%% Define the nodes of the beam
Node1=en_free([0,0,0]);      % Where the node begins
Node2=en_free([6 0 0]);     % Where the section from conical becomes cylindrical
Node3=en_free([18 0 0]);    % Where the wings are clamped
Node4=en_free([40 0 0]);    % Where the conical section of the tail begins
Node5=en_free([48 0 2.5]);  % Where the rudder is clamped
Node6=en_free([52 0 2.5]);  % Where the aircraft ends

aircraft.en=[Node1, Node2, Node3, Node4, Node5, Node6];

% clear unuseful variables
clear Node1
clear Node2
clear Node3
clear Node4
clear Node5
clear Node6
clear i

%% Discretisation parameters
nel_nose=18*2;
nel_front_fuselage=35;
nel_rear_fuselage=63;
nel_front_tail=23*2;
nel_rear_tail=11*2;

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
