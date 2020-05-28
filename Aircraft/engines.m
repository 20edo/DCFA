% This script builds the beam that represent the engines supports and adds
% the properties of the engines
% The origin is set at the nose of the aircraft, the x axis is align with
% the fuselage ond points to the front

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

%% Discretisation parameters
nel=10;         % Number of elements

%% Material properties (Alluminium)

E=70*1e9;        % Young modulus
G=27*1e9;        % Shear modulus
rho=2700;    % Density

%% Throttle

t1=1;           % Percentage of thrust engine 1 wrt maximum at sea level
t2=1;           % Percentage of thrust engine 2 wrt maximum at sea level
t3=1;           % Percentage of thrust engine 3 wrt maximum at sea level
t4=1;           % Percentage of thrust engine 4 wrt maximum at sea level

%% Engine properties and positions

T0 = 179.9; %KN % Nominal thrust
Me = 3220; %Kg  % Mass engine
De = 2.15; %m % External diameter
Re = 2.15/2; %m % External radius
Le = 3.73; %m
xp = Re/2; %m
yp = Re/2; %m
zp = Le/2; %m
Ix = (Me*Re^2)/2 + Me*(yp^2 + zp^2); %m^4
Iy = Me/12*(3*Re^2 + Le^2) + Me*(xp^2 + zp^2); %m^4
Iz = Me/12*(3*Re^2 + Le^2) + Me*(xp^2 + yp^2); %m^4

M = [Me 0 0 0 0 0;
    0 Me 0 0 0 0;
    0 0 Me 0 0 0;
    0 0 0 Ix 0 0;
    0 0 0 0 Iy 0;
    0 0 0 0 0 Iz];
K = T0*[0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 1 0 0; % dw*T0*theta
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

x=[-5.5 0 -2]';                      % Relative position of the engine wrt nodes in the wings

Node17=en_free(aircraft.en(8).x+x);
Node18=en_free(aircraft.en(9).x+x);
Node19=en_free(aircraft.en(11).x+x);
Node20=en_free(aircraft.en(12).x+x);

Engines=[Node17 Node18 Node19 Node20];

for i=1:length(Engines)
    Engines(i).K=zeros(6);          % K matrix of the engine
    Engines(i).M=zeros(6);          % Inertia matrix of the engine
    aircraft.en=[aircraft.en Engines(i)];
end

%% Supports

% Support properties

c=@(x) 1+0.*x;  
h=@(x) 0.12+0.*x; 
t=@(x) 0.020+0.*x;

% support 1

L=norm(aircraft.en(17).x-aircraft.en(8).x,2);
support_1=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_1.o=aircraft.en(8).x;
support_1.vx=aircraft.en(17).x-aircraft.en(8).x;
support_1.vx=support_1.vx/norm(support_1.vx,2);
support_1.vy=[0 0 1]';
support_1.vy=support_1.vy-dot(support_1.vy,support_1.vx)*support_1.vx;
support_1.vy=support_1.vy/norm(support_1.vy);
support_1.name='support_1';

% support 2

L=norm(aircraft.en(18).x-aircraft.en(9).x,2);
support_2=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_2.o=aircraft.en(9).x;
support_2.vx=aircraft.en(18).x-aircraft.en(9).x;
support_2.vx=support_2.vx/norm(support_2.vx,2);
support_2.vy=[0 0 1]';
support_2.vy=support_2.vy-dot(support_2.vy,support_2.vx)*support_2.vx;
support_2.vy=support_2.vy/norm(support_2.vy);
support_2.name='support_2';

% support 3

L=norm(aircraft.en(19).x-aircraft.en(11).x,2);
support_3=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_3.o=aircraft.en(11).x;
support_3.vx=aircraft.en(19).x-aircraft.en(11).x;
support_3.vx=support_3.vx/norm(support_3.vx,2);
support_3.vy=[0 0 1]';
support_3.vy=support_3.vy-dot(support_3.vy,support_3.vx)*support_3.vx;
support_3.vy=support_3.vy/norm(support_3.vy);
support_3.name='support_3';

% support 4

L=norm(aircraft.en(20).x-aircraft.en(12).x,2);
support_4=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_4.o=aircraft.en(12).x;
support_4.vx=aircraft.en(20).x-aircraft.en(12).x;
support_4.vx=support_4.vx/norm(support_4.vx,2);
support_4.vy=[0 0 1]';
support_4.vy=support_4.vy-dot(support_4.vy,support_4.vx)*support_4.vx;
support_4.vy=support_4.vy/norm(support_4.vy);
support_4.name='support_4';


%% Add supports to the model

support=[ support_1 support_2 support_3 support_4];

for i=1:length(support)
    aircraft=m_add_beam(aircraft,support(i));
end

%% Clear unuseful variables
clear c
clear De
clear Dx
clear E
clear Engines
clear G
clear i
clear h
clear Ix
clear Iy
clear Iz 
clear K
clear L
clear Le
clear M
clear Me
clear nel
clear Node17
clear Node18
clear Node19
clear Node20
clear Re
clear rho
clear support_1
clear support_2
clear support_3
clear support_4
clear t
clear T0
clear t1
clear t2
clear t3
clear t4
clear x
clear xp
clear yp
clear zp
clear support