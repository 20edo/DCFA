% This script builds the beam that represent the wings of the aircraft
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
nel_1=30;
nel_2=30;
nel_3=30;


%% Parameters of the wing

M=eye(6);                                       % Mass matrix of the engine
x1=[0 10 0]';                                   % Relative position of the engine 1 (right) wrt node 3
x2=[0 20 0]';                                   % Relative position of the engine 2 (right) wrt node 3
x3=[0 30 0]';                                   % Relative position of the tip of the wing (right) wrt node 3

sq=[1.2 1 0.8];                                 % latus of the squares of the different sections of the wing
Node7=en_mass(aircraft.en(3).x+x1,M);           % Position of the right engine closer to the root
Node8=en_mass(aircraft.en(3).x+x2,M);           % Position of the right engine farther from the root
Node9=en_free(aircraft.en(3).x+x3);             % Position of the tip of the right wing

Node10=en_mass(aircraft.en(3).x+x1-[0 2*x1(2) 0]',M);           % Position of the left engine closer to the root
Node11=en_mass(aircraft.en(3).x+x2-[0 2*x2(2) 0]',M);           % Position of the left engine farther from the root
Node12=en_free(aircraft.en(3).x+x3-[0 2*x3(2) 0]');             % Position of the tip of the left wing

aircraft.en=[aircraft.en Node7 Node8 Node9 Node10 Node11 Node12];

%% Material properties (Aluminium)

E=70*1e3;   % Young modulus
G=27*1e3;   % Shear modulus
rho=2700;   % Density


%% Right wing

% First beam of the right wing

L=norm(aircraft.en(7).x-aircraft.en(3).x,2);
l=sq(1);
r_wing1=b_constant_p_square(L,l,E,G,rho,nel_1);
r_wing1.o=aircraft.en(3).x;
r_wing1.vx=aircraft.en(7).x-aircraft.en(3).x;
r_wing1.vx=r_wing1.vx/norm(r_wing1.vx,2);
r_wing1.vy=[0 1 0]';
r_wing1.name='r_wing1';

% Second beam of the right wing

L=norm(aircraft.en(8).x-aircraft.en(7).x,2);
l=sq(2);
r_wing2=b_constant_p_square(L,l,E,G,rho,nel_2);
r_wing2.o=aircraft.en(7).x;
r_wing2.vx=aircraft.en(8).x-aircraft.en(7).x;
r_wing2.vx=r_wing2.vx/norm(r_wing2.vx,2);
r_wing2.vy=[0 1 0]';
r_wing2.name='r_wing2';


% Third beam of the right wing

L=norm(aircraft.en(9).x-aircraft.en(8).x,2);
l=sq(3);
r_wing3=b_constant_p_square(L,l,E,G,rho,nel_3);
r_wing3.o=aircraft.en(8).x;
r_wing3.vx=aircraft.en(9).x-aircraft.en(8).x;
r_wing3.vx=r_wing3.vx/norm(r_wing3.vx,2);
r_wing3.vy=[0 1 0]';
r_wing3.name='r_wing3';

%% Left wing 

% First beam of the left wing

L=norm(aircraft.en(10).x-aircraft.en(3).x,2);
l=sq(1);
l_wing1=b_constant_p_square(L,l,E,G,rho,nel_1);
l_wing1.o=aircraft.en(3).x;
l_wing1.vx=aircraft.en(10).x-aircraft.en(3).x;
l_wing1.vx=l_wing1.vx/norm(l_wing1.vx,2);
l_wing1.vy=[0 1 0]';
l_wing1.name='l_wing1';

% Second beam of the left wing

L=norm(aircraft.en(11).x-aircraft.en(10).x,2);
l=sq(2);
l_wing2=b_constant_p_square(L,l,E,G,rho,nel_2);
l_wing2.o=aircraft.en(10).x;
l_wing2.vx=aircraft.en(11).x-aircraft.en(10).x;
l_wing2.vx=l_wing2.vx/norm(l_wing2.vx,2);
l_wing2.vy=[0 1 0]';
l_wing2.name='l_wing2';

% Third beam of the left wing

L=norm(aircraft.en(12).x-aircraft.en(11).x,2);
l=sq(3);
l_wing3=b_constant_p_square(L,l,E,G,rho,nel_3);
l_wing3.o=aircraft.en(11).x;
l_wing3.vx=aircraft.en(12).x-aircraft.en(11).x;
l_wing3.vx=l_wing3.vx/norm(l_wing3.vx,2);
l_wing3.vy=[0 1 0]';
l_wing3.name='l_wing3';

%% Add beams to the aircraft

wings_beams=[ r_wing1 r_wing2 r_wing3 l_wing1 l_wing2 l_wing3];

for i=1:length(wings_beams)
    aircraft=m_add_beam(aircraft,wings_beams(i));
end

%% Clear unusefull variables

clear sq
clear E
clear G
clear i
clear l
clear L
clear l_wing1
clear l_wing2
clear l_wing3
clear M
clear nel_1
clear nel_2
clear nel_3
clear Node10
clear Node11
clear Node12
clear Node7
clear Node8
clear Node9
clear r_wing1
clear r_wing2
clear r_wing3
clear rho
clear x1
clear x2
clear x3



