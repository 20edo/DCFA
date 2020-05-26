
cd ..
cd ..

addpath('Section')
addpath('Element')
addpath('Beam')
addpath('internal_nodes')
addpath('External_nodes')
addpath('Model')

cd Aircraft\


%% Initialize model

aircraft=m_init();

Node1=en_ground([0,0,0]);      % Where the node begins

aircraft.en=[Node1];

%% Discretisation parameters
nel_1=30;
nel_2=30;
nel_3=30;


%% Parameters of the wing

x0=[0 0 2]';                                    % Higth of the wings
x1=[-4.6 10 0]';                                % Relative position of the engine 1 (right) wrt node 3
x1(3)=-norm(x1)*tand(3);                         % Added anhedral angle
x2=[-7  15 0]';                                 % Relative position of the engine 2 (right) wrt node 3
x2(3)=-norm(x2)*tand(3);                         % Added anhedral angle
x3=[-13 27 0]';                                 % Relative position of the tip of the wing (right) wrt node 3
x3(3)=-norm(x3)*tand(3);                         % Added anhedral angle
i_ang=2;                                        % Angle of incidence (calettamento)


chord=@(x) 7.72+2-x/norm(x3)*5;            % Chord of the wing
heigth=@(x) 0.12+0.*x;                          % Relative heigth of the wing section
thickness=@(x) 0.020+0.01-x/norm(x3)*0.02;      % Thickness of the panels (see ssh_profile) 

Node2=en_free(Node1.x+x1);                      % Position of the right engine closer to the root
Node3=en_free(Node1.x+x2);                      % Position of the right engine farther from the root
Node4=en_free(Node1.x+x3);                     % Position of the tip of the right wing

Node5=en_free(Node1.x+x1-[0 2*x1(2) 0]');      % Position of the left engine closer to the root
Node6=en_free(Node1.x+x2-[0 2*x2(2) 0]');      % Position of the left engine farther from the root
Node7=en_free(Node1.x+x3-[0 2*x3(2) 0]');      % Position of the tip of the left wing

aircraft.en=[aircraft.en Node2 Node3 Node4 Node5 Node6 Node7];

%% Material properties (Aluminium)

E=70*1e9;   % Young modulus
G=27*1e9;   % Shear modulus
rho=2700;   % Density


%% Right wing

% First beam of the left wing

L=norm(aircraft.en(2).x-aircraft.en(1).x,2);
c=@(x) chord(x);
t=@(x) thickness(x);
h=@(x) heigth(x);
l_wing1=b_ssh_profile(L,c,h,t,E,G,rho,nel_1);
l_wing1.o=aircraft.en(1).x;
l_wing1.vx=aircraft.en(2).x-aircraft.en(1).x;
l_wing1.vx=l_wing1.vx/norm(l_wing1.vx,2);
l_wing1.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
l_wing1.vy=l_wing1.vy-dot(l_wing1.vy,l_wing1.vx)*l_wing1.vx;
l_wing1.vy=l_wing1.vy/norm(l_wing1.vy);
l_wing1.name='l_wing1';

% Second beam of the left wing

L=norm(aircraft.en(3).x-aircraft.en(2).x,2);
c=@(x) chord(x+norm(x1));
t=@(x) thickness(x+norm(x1));
h=@(x) heigth(x+norm(x1));
l_wing2=b_ssh_profile(L,c,h,t,E,G,rho,nel_2);
l_wing2.o=aircraft.en(2).x;
l_wing2.vx=aircraft.en(3).x-aircraft.en(2).x;
l_wing2.vx=l_wing2.vx/norm(l_wing2.vx,2);
l_wing2.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
l_wing2.vy=l_wing2.vy-dot(l_wing2.vy,l_wing2.vx)*l_wing2.vx;
l_wing2.vy=l_wing2.vy/norm(l_wing2.vy);
l_wing2.name='l_wing2';


% Third beam of the left wing


L=norm(aircraft.en(4).x-aircraft.en(3).x,2);
c=@(x) chord(x+norm(x2));
t=@(x) thickness(x+norm(x2));
h=@(x) heigth(x+norm(x2));
l_wing3=b_ssh_profile(L,c,h,t,E,G,rho,nel_3);
l_wing3.o=aircraft.en(3).x;
l_wing3.vx=aircraft.en(4).x-aircraft.en(3).x;
l_wing3.vx=l_wing3.vx/norm(l_wing3.vx,2);
l_wing3.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
l_wing3.vy=l_wing3.vy-dot(l_wing3.vy,l_wing3.vx)*l_wing3.vx;
l_wing3.vy=l_wing3.vy/norm(l_wing3.vy);
l_wing3.name='l_wing3';

%% Left wing 

% First beam of the right wing

L=norm(aircraft.en(5).x-aircraft.en(1).x,2);
c=@(x) chord(x);
t=@(x) thickness(x);
h=@(x) heigth(x);
r_wing1=b_ssh_profile(L,c,h,t,E,G,rho,nel_1);
r_wing1.o=aircraft.en(1).x;
r_wing1.vx=aircraft.en(5).x-aircraft.en(1).x;
r_wing1.vx=r_wing1.vx/norm(r_wing1.vx);
r_wing1.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
r_wing1.vy=r_wing1.vy-dot(r_wing1.vy,l_wing1.vx)*r_wing1.vx;
r_wing1.vy=r_wing1.vy/norm(r_wing1.vy);
r_wing1.name='r_wing1';

% Second beam of the right wing

L=norm(aircraft.en(6).x-aircraft.en(5).x,2);
c=@(x) chord(x+norm(x1));
t=@(x) thickness(x+norm(x1));
h=@(x) heigth(x+norm(x1));
r_wing2=b_ssh_profile(L,c,h,t,E,G,rho,nel_2);
r_wing2.o=aircraft.en(5).x;
r_wing2.vx=aircraft.en(6).x-aircraft.en(5).x;
r_wing2.vx=r_wing2.vx/norm(r_wing2.vx,2);
r_wing2.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
r_wing2.vy=r_wing2.vy-dot(r_wing2.vy,r_wing2.vx)*r_wing2.vx;
r_wing2.vy=r_wing2.vy/norm(r_wing2.vy);
r_wing2.name='r_wing2';

% Third beam of the left wing

L=norm(aircraft.en(7).x-aircraft.en(6).x,2);
c=@(x) chord(x+norm(x2));
t=@(x) thickness(x+norm(x2));
h=@(x) heigth(x+norm(x2));
r_wing3=b_ssh_profile(L,c,h,t,E,G,rho,nel_3);
r_wing3.o=aircraft.en(6).x;
r_wing3.vx=aircraft.en(7).x-aircraft.en(6).x;
r_wing3.vx=r_wing3.vx/norm(r_wing3.vx,2);
r_wing3.vy=[-cosd(i_ang) 0 -sind(i_ang)]';
r_wing3.vy=r_wing3.vy-dot(r_wing3.vy,r_wing3.vx)*r_wing3.vx;
r_wing3.vy=r_wing3.vy/norm(r_wing3.vy);
r_wing3.name='r_wing3';

%% Add beams to the aircraft

wings_beams=[l_wing1 l_wing2 l_wing3 r_wing1 r_wing2 r_wing3];

for i=1:length(wings_beams)
    aircraft=m_add_beam(aircraft,wings_beams(i));
%     disp(wings_beams(i).name)
%     M=sum(sum(wings_beams(i).M(1:6:end,1:6:end));
%     disp('M=')
%     disp(M)
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
clear wings_beams
clear x0
clear wing_support
clear t
clear thickness
clear h
clear c
clear chord
clear heigth
clear Node13
clear i_ang

aircraft=m_compute_matrices(aircraft);

if 0
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    [fig] = m_plot3d(aircraft,options);
end

cd .. 
cd Test\Playground