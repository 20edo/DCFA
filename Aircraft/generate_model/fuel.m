% This script edits the beam that represent the wings of the aircraft to
% include the fuel tanks
% The origin is set at the nose of the aircraft, the x axis is aligned with
% the fuselage and points to the front
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

%% Material properties (JET fuel A1)

E=0;        % Young modulus
G=0;        % Shear modulus
rho=804;    % Density

%% Fuel percentage
study = 4; 
% 1 -> reference condition 
% 2 -> 100% of fuel 
% 3 -> no fuel 
% 4 -> fuel only in the internal tanks
% 5 -> fuel only in the external tanks 
%
switch study
    case 1 
% Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=1;      % Percentage of fuel tank1
fuel_t2=1;      % Percentage of fuel tank2
fuel_t3=0.5;      % Percentage of fuel tank3
% Left
fuel_t4=1;      % Percentage of fuel tank4
fuel_t5=1;      % Percentage of fuel tank5
fuel_t6=0.5;      % Percentage of fuel tank6
    case 2 
       % Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=1;      % Percentage of fuel tank1
fuel_t2=1;      % Percentage of fuel tank2
fuel_t3=1;      % Percentage of fuel tank3
% Left
fuel_t4=1;      % Percentage of fuel tank4
fuel_t5=1;      % Percentage of fuel tank5
fuel_t6=1;      % Percentage of fuel tank6 
    case 3 
        % Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=0;      % Percentage of fuel tank1
fuel_t2=0;      % Percentage of fuel tank2
fuel_t3=0;      % Percentage of fuel tank3
% Left
fuel_t4=0;      % Percentage of fuel tank4
fuel_t5=0;      % Percentage of fuel tank5
fuel_t6=0;      % Percentage of fuel tank6
    case 4 
        % Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=1;      % Percentage of fuel tank1
fuel_t2=1;      % Percentage of fuel tank2
fuel_t3=0;      % Percentage of fuel tank3
% Left
fuel_t4=1;      % Percentage of fuel tank4
fuel_t5=1;      % Percentage of fuel tank5
fuel_t6=0;      % Percentage of fuel tank6
    case 5 
% Tanks are ordered from the fuselage to the tip of the wing
% Right
fuel_t1=1;      % Percentage of fuel tank1
fuel_t2=1;      % Percentage of fuel tank2
fuel_t3=0.3;      % Percentage of fuel tank3
% Left
fuel_t4=1;      % Percentage of fuel tank4
fuel_t5=1;      % Percentage of fuel tank5
fuel_t6=0.3;      % Percentage of fuel tank6
end
% %% Study of the change of the fuel 
% load carb.mat
% i = carb(end,1); 
% fuel_t1=carb(i,1);      % Percentage of fuel tank1
% fuel_t2=carb(i,2);      % Percentage of fuel tank2
% fuel_t3=carb(i,3);      % Percentage of fuel tank3
% % Left
% fuel_t4=carb(i,1);      % Percentage of fuel tank4
% fuel_t5=carb(i,2);      % Percentage of fuel tank5
% fuel_t6=carb(i,3);      % Percentage of fuel tank6
% i = i+1; 
% carb(end,1) = i; 
% save carb


%% Tank
% Tank1

L=norm(aircraft.b(7).L);
l1=1.6;                        % Dimension of the tank along the chord
l2=0.16;                       % Heigth of the tank
nel=(size(aircraft.b(7).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(7).o;
tank.vx=aircraft.b(7).vx;
tank.vy=aircraft.b(7).vy;
% disp('Tank_1')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(7).M=aircraft.b(7).M+fuel_t1*tank.M;

% Tank2

L=norm(aircraft.b(8).L);
l1=1.4;
l2=0.14;
nel=(size(aircraft.b(8).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(8).o;
tank.vx=aircraft.b(8).vx;
tank.vy=aircraft.b(8).vy;
% disp('Tank_2')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(8).M=aircraft.b(8).M+fuel_t2*tank.M;

% Tank3

L=norm(aircraft.b(9).L);
l1=1.0;
l2=0.10;
nel=(size(aircraft.b(9).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(9).o;
tank.vx=aircraft.b(9).vx;
tank.vy=aircraft.b(9).vy;
% disp('Tank_3')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(9).M=aircraft.b(9).M+fuel_t3*tank.M;

% Tank4

L=norm(aircraft.b(10).L);
l1=1.6;                        % Dimension of the tank along the chord
l2=0.16;                       % Heigth of the tank
nel=(size(aircraft.b(10).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(10).o;
tank.vx=aircraft.b(10).vx;
tank.vy=aircraft.b(10).vy;
% disp('Tank_4')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(10).M=aircraft.b(10).M+fuel_t4*tank.M;

% Tank5

L=norm(aircraft.b(11).L);
l1=1.4;
l2=0.14;
nel=(size(aircraft.b(11).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(11).o;
tank.vx=aircraft.b(11).vx;
tank.vy=aircraft.b(11).vy;
% disp('Tank_5')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(11).M=aircraft.b(11).M+fuel_t5*tank.M;

% Tank6

L=norm(aircraft.b(12).L);
l1=1;
l2=0.10;
nel=(size(aircraft.b(12).M,1)-6)/6;
tank=b_constant_p_rect(L,l1,l2,E,G,rho,nel);
tank.o=aircraft.b(12).o;
tank.vx=aircraft.b(12).vx;
tank.vy=aircraft.b(12).vy;
% disp('Tank_6')
% disp(sum(sum(tank.M(1:6:end,1:6:end))))
aircraft.b(12).M=aircraft.b(12).M+fuel_t6*tank.M;




%% Clear unuseful variables

clear fuel_t1
clear fuel_t2
clear fuel_t3
clear fuel_t4
clear fuel_t5
clear fuel_t6
clear l1
clear l2
clear nel
clear tank
clear E
clear G
clear L
clear rho

