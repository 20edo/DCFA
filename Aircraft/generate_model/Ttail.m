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

E=70*1e9;   % Young modulus
G=27*1e9;   % Shear modulus
rho=2700;   % Density

%% Parameters of the tail

x1=[6 0 9.6]';     % Position of the stabilizer wrt node 5
x2=[4 8.7 0]';     % Position of the tip of the stabilizer wrt the center
x2(3)=-norm(x1)*tand(3);                         % Added anhedral angle


i_ang=-5;                                        % Angle of incidence (calettamento)

Node14=en_free(aircraft.en(5).x+x1);
Node15=en_free(aircraft.en(5).x+x1+x2);
Node16=en_free(aircraft.en(5).x+x1+x2-[0 2*x2(2) 0]');

aircraft.en=[aircraft.en Node14 Node15 Node16];

%% Rudder

L=norm(aircraft.en(14).x-aircraft.en(5).x,2);
c=@(x) 5+0.*x;  
h=@(x) 0.12+0.*x; 
t=@(x) 0.020+0.*x;
rudder=b_ssh_profile(L,c,h,t,E,G,rho,nel_1);
rudder.o=aircraft.en(5).x;
rudder.vx=aircraft.en(14).x-aircraft.en(5).x;
rudder.vx=rudder.vx/norm(rudder.vx,2);
rudder.vy=[-1 0 0]';
rudder.vy=rudder.vy-dot(rudder.vy,rudder.vx)*rudder.vx;
rudder.vy=rudder.vy/norm(rudder.vy);
rudder.name='rudder';



%% Stabilizer

% Right


L=norm(aircraft.en(15).x-aircraft.en(14).x,2);
c=@(x) 4-2.*x/L;  
h=@(x) 0.12+0.*x; 
t=@(x) 0.010+0.*x;
r_stabilizer=b_ssh_profile(L,c,h,t,E,G,rho,nel_2);
r_stabilizer.o=aircraft.en(14).x;
r_stabilizer.vx=aircraft.en(15).x-aircraft.en(14).x;
r_stabilizer.vx=r_stabilizer.vx/norm(r_stabilizer.vx,2);
r_stabilizer.vy=[cosd(i_ang) 0 -sind(i_ang)]';
r_stabilizer.vy=r_stabilizer.vy-dot(r_stabilizer.vy,r_stabilizer.vx)*r_stabilizer.vx;
r_stabilizer.vy=r_stabilizer.vy/norm(r_stabilizer.vy);
r_stabilizer.name='r_stabilizer';

% Left

L=norm(aircraft.en(16).x-aircraft.en(14).x,2);
c=@(x) 4-2*x/L; 
h=@(x) 0.12+0.*x; 
t=@(x) 0.010+0.*x;
l_stabilizer=b_ssh_profile(L,c,h,t,E,G,rho,nel_3);
l_stabilizer.o=aircraft.en(14).x;
l_stabilizer.vx=aircraft.en(16).x-aircraft.en(14).x;
l_stabilizer.vx=l_stabilizer.vx/norm(l_stabilizer.vx,2);
l_stabilizer.vy=[cosd(i_ang) 0 -sind(i_ang)]';
l_stabilizer.vy=l_stabilizer.vy-dot(l_stabilizer.vy,l_stabilizer.vx)*l_stabilizer.vx;
l_stabilizer.vy=l_stabilizer.vy/norm(l_stabilizer.vy);
l_stabilizer.name='r_stabilizer';

%% Add beams to the model

Ttail_beams=[rudder r_stabilizer l_stabilizer];

for i=1:length(Ttail_beams)
    aircraft=m_add_beam(aircraft,Ttail_beams(i));
%     disp(Ttail_beams(i).name)
%     M=sum(sum(Ttail_beams(i).M(1:6:end,1:6:end)));
%     disp('M=')
%     disp(M)

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
clear Node16
clear Node14
clear Node15
clear r_stabilizer
clear rho
clear rudder
clear x1
clear x2
clear Ttail_beams
clear c 
clear h
clear t
clear i_ang

