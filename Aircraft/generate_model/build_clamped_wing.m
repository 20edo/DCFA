function model=build_clamped_wing(aircraft,nel)
% This function builds a model that represent the wing clamped.
% aircraft  ->      model of the aircraft
% nel       ->      number of elements

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

engine=true;

%% Discretization parameter
Ltot=aircraft.b(7).L+aircraft.b(8).L+aircraft.b(9).L;
nel_1=ceil(nel*aircraft.b(7).L/Ltot);
nel_2=ceil(nel*aircraft.b(8).L/Ltot);
nel_3=nel-nel_1-nel_2;

%% Wing properties

E=70*1e9;   % Young modulus
G=27*1e9;   % Shear modulus
rho=2700;   % Density


%% Init model and add ground node

model=m_init();
if engine
    model.en=[en_ground(aircraft.en(7).x) ...
        en_mass(aircraft.en(17).x,aircraft.en(17).M) en_mass(aircraft.en(18).x,aircraft.en(18).M)];
else
    model.en=[en_ground(aircraft.en(7).x)];
end
%% First beam

L=aircraft.b(7).L;
c=@(x) aircraft.b(7).c(x);
t=@(x) aircraft.b(7).t(x);
h=@(x) aircraft.b(7).h(x);
r_wing1=b_ssh_profile(L,c,h,t,E,G,rho,nel_1);
r_wing1.o=aircraft.b(7).o;
r_wing1.vx=aircraft.b(7).vx;
r_wing1.vy=aircraft.b(7).vy;
r_wing1.name='r_wing1';


%% Second beam

L=aircraft.b(8).L;
c=@(x) aircraft.b(8).c(x);
t=@(x) aircraft.b(8).t(x);
h=@(x) aircraft.b(8).h(x);
r_wing2=b_ssh_profile(L,c,h,t,E,G,rho,nel_2);
r_wing2.o=aircraft.b(8).o;
r_wing2.vx=aircraft.b(8).vx;
r_wing2.vy=aircraft.b(8).vy;
r_wing2.name='r_wing2';

%% Third beam

L=aircraft.b(9).L;
c=@(x) aircraft.b(9).c(x);
t=@(x) aircraft.b(9).t(x);
h=@(x) aircraft.b(9).h(x);
r_wing3=b_ssh_profile(L,c,h,t,E,G,rho,nel_3);
r_wing3.o=aircraft.b(9).o;
r_wing3.vx=aircraft.b(9).vx;
r_wing3.vy=aircraft.b(9).vy;
r_wing3.name='r_wing3';

%% Add fuel

% Material properties (JET fuel A1)

E=0;        % Young modulus
G=0;        % Shear modulus
rho=804;    % Density

% Fuel percentage
% Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=1;      % Percentage of fuel tank1
fuel_t2=1;      % Percentage of fuel tank2
fuel_t3=1;      % Percentage of fuel tank3
% Left
fuel_t4=1;      % Percentage of fuel tank4
fuel_t5=1;      % Percentage of fuel tank5
fuel_t6=1;      % Percentage of fuel tank6

% Tank
% Tank1

L=norm(aircraft.b(7).L);
l1=1.6;                        % Dimension of the tank along the chord
l2=0.16;                       % Heigth of the tank
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel_1);
tank.o=aircraft.b(7).o;
tank.vx=aircraft.b(7).vx;
tank.vy=aircraft.b(7).vy;
% disp('Tank_1')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
r_wing1.M=r_wing1.M+fuel_t1*tank.M;

% Tank2

L=norm(aircraft.b(8).L);
l1=1.4;
l2=0.14;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel_2);
tank.o=aircraft.b(8).o;
tank.vx=aircraft.b(8).vx;
tank.vy=aircraft.b(8).vy;
% disp('Tank_2')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
r_wing2.M=r_wing2.M+fuel_t2*tank.M;

% Tank3

L=norm(aircraft.b(9).L);
l1=1.0;
l2=0.10;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel_3);
tank.o=aircraft.b(9).o;
tank.vx=aircraft.b(9).vx;
tank.vy=aircraft.b(9).vy;
% disp('Tank_3')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
r_wing3.M=r_wing3.M+fuel_t3*tank.M;




%% Add beams to model

if engine
    % The last two are the engine support
    wing=[r_wing1 r_wing2 r_wing3 aircraft.b(16) aircraft.b(17)];
else
    wing=[r_wing1 r_wing2 r_wing3];
end
    
for i=1:length(wing)
    model=m_add_beam(model,wing(i));
end

