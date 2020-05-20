% This script builds the T tail of the aircraft

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

%% Material properties (Aluminium)

E=70*1e3;   % Young modulus
G=27*1e3;   % Shear modulus
rho=2700;   % Density

%% Parameters of the tail

x1=[-0.5 0 1]';      % Position of the stabilizer wrt node 5
x2=[-1 1 0]';        % Position of the tip of the stabilizer wrt the center
l=0.3;               % Latus of the square

Node13=en_free(aircraft.en(5).x+x1);
Node14=en_free(aircraft.en(5).x+x1+x2);
Node15=en_free(aircraft.en(5).x+x1+x2-[0 2*x2(2) 0]');

aircraft.en=[aircraft.en Node13 Node14 Node15];

%% Rudder

L=norm(aircraft.en(13).x-aircraft.en(5).x,2);
rudder=b_constant_p_square(L,l,E,G,rho,nel_3);
rudder.o=aircraft.en(5).x;
rudder.vx=aircraft.en(13).x-aircraft.en(5).x;
rudder.vx=rudder.vx/norm(rudder.vx,2);
rudder.vy=[0 1 0]';
rudder.name='rudder';



%% Stabilizer

% Right 


L=norm(aircraft.en(14).x-aircraft.en(13).x,2);
r_stabilizer=b_constant_p_square(L,l,E,G,rho,nel_3);
r_stabilizer.o=aircraft.en(13).x;
r_stabilizer.vx=aircraft.en(14).x-aircraft.en(13).x;
r_stabilizer.vx=r_stabilizer.vx/norm(r_stabilizer.vx,2);
r_stabilizer.vy=[0 1 0]';
r_stabilizer.name='r_stabilizer';

% Left

L=norm(aircraft.en(15).x-aircraft.en(13).x,2);
l_stabilizer=b_constant_p_square(L,l,E,G,rho,nel_3);
l_stabilizer.o=aircraft.en(13).x;
l_stabilizer.vx=aircraft.en(15).x-aircraft.en(13).x;
l_stabilizer.vx=l_stabilizer.vx/norm(l_stabilizer.vx,2);
l_stabilizer.vy=[0 1 0]';
l_stabilizer.name='l_stabilizer';

%% Add beams to the model

Ttail_beams=[rudder r_stabilizer l_stabilizer];

for i=1:length(Ttail_beams)
    aircraft=m_add_beam(aircraft,Ttail_beams(i));
end

clear E
clear G
clear i
clear l
clear L
clear l_stabilizer
clear nel_1
clear nel_2
clear nel_3
clear Node13
clear Node14
clear Node15
clear r_stabilizer
clear rho
clear rudder
clear x1
clear x2
