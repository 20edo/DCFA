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

E=70*1e3;        % Young modulus
G=27*1e3;        % Shear modulus
rho=2700;    % Density

%% Throttle

t1=1;           % Percentage of thrust engine 1 wrt maximum at sea level
t2=1;           % Percentage of thrust engine 2 wrt maximum at sea level
t3=1;           % Percentage of thrust engine 3 wrt maximum at sea level
t4=1;           % Percentage of thrust engine 4 wrt maximum at sea level

%% Engine properties and positions

x=[5.5 0 -2]';                      % Relative position of the engine wrt nodes in the wings

Node18=en_free(aircraft.en(8).x+x);
Node19=en_free(aircraft.en(9).x+x);
Node20=en_free(aircraft.en(11).x+x);
Node21=en_free(aircraft.en(12).x+x);

Engines=[Node18 Node19 Node20 Node21];

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

L=norm(aircraft.en(14).x-aircraft.en(5).x,2);
support_1=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_1.o=aircraft.en(5).x;
support_1.vx=aircraft.en(18).x-aircraft.en(8).x;
support_1.vx=support_1.vx/norm(support_1.vx,2);
support_1.vy=[0 0 1]';
support_1.name='support_1';

% support 2

L=norm(aircraft.en(19).x-aircraft.en(9).x,2);
support_2=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_2.o=aircraft.en(9).x;
support_2.vx=aircraft.en(19).x-aircraft.en(9).x;
support_2.vx=support_2.vx/norm(support_2.vx,2);
support_2.vy=[0 0 1]';
support_2.name='support_2';

% support 3

L=norm(aircraft.en(20).x-aircraft.en(11).x,2);
support_3=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_3.o=aircraft.en(11).x;
support_3.vx=aircraft.en(20).x-aircraft.en(11).x;
support_3.vx=support_3.vx/norm(support_3.vx,2);
support_3.vy=[0 0 1]';
support_3.name='support_3';

% support 4

L=norm(aircraft.en(21).x-aircraft.en(12).x,2);
support_4=b_ssh_profile(L,c,h,t,E,G,rho,nel);
support_4.o=aircraft.en(12).x;
support_4.vx=aircraft.en(21).x-aircraft.en(12).x;
support_4.vx=support_4.vx/norm(support_4.vx,2);
support_4.vy=[0 0 1]';
support_4.name='support_4';