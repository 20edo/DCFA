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
model.en=[en_ground(aircraft.en(7).x) aircraft.en(17) aircraft.en(18)];

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
r_wing3=b_ssh_profile(L,c,h,t,E,G,rho,nel_1);
r_wing3.o=aircraft.b(9).o;
r_wing3.vx=aircraft.b(9).vx;
r_wing3.vy=aircraft.b(9).vy;
r_wing3.name='r_wing3';

%% Add beams to model

% The last two are the engine support
wing=[r_wing1 r_wing2 r_wing3 aircraft.b(16) aircraft.b(17)];

for i=1:length(wing)
    model=m_add_beam(model,wing(i));
end
