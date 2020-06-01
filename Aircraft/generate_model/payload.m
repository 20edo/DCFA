% This script edits the beam that represent the fuselage of the aircraft to
% include the payload
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

%% Kind of payload
human=true;
vehicles=false;

%% Material properties (Instruments)

E=0;            % Young modulus
G=0;            % Shear modulus
rho=2700/20;    % Density

%% Payload

% Cockpit and instruments

L=norm(aircraft.b(1).L);
l=1;
nel=(size(aircraft.b(1).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(1).o;
Payload.vx=aircraft.b(1).vx;
Payload.vy=aircraft.b(1).vy;
% disp('Cockpit')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(1).M=aircraft.b(1).M+Payload.M;

% Nose instuments

L=norm(aircraft.b(2).L);
l=1;
nel=(size(aircraft.b(2).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(2).o;
Payload.vx=aircraft.b(2).vx;
Payload.vy=aircraft.b(2).vy;
% disp('Nose ')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(2).M=aircraft.b(2).M+Payload.M;

if human
    
    % Material properties (Soldiers)

    E=0;            % Young modulus
    G=0;            % Shear modulus
    rho=1000/10;    % Density


    % Soldiers    
    L=norm(aircraft.b(3).L);
    l=2;
    nel=(size(aircraft.b(3).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(3).o;
    Payload.vx=aircraft.b(3).vx;
    Payload.vy=aircraft.b(3).vy;
%     disp('Soldiers_1')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(3).M=aircraft.b(3).M+Payload.M;
    
    
    % Soldiers    
    L=norm(aircraft.b(4).L);
    l=2;
    nel=(size(aircraft.b(4).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(4).o;
    Payload.vx=aircraft.b(4).vx;
    Payload.vy=aircraft.b(4).vy;
%     disp('Soldiers_2')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(4).M=aircraft.b(4).M+Payload.M;
    
    
elseif vehicles
    
    %% Material properties (Vehicles)

    E=0;            % Young modulus
    G=0;            % Shear modulus
    rho=8000/20;    % Density
    
    % Vehicles
    
    L=norm(aircraft.b(3).L);
    l=2.5;
    nel=(size(aircraft.b(3).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(3).o;
    Payload.vx=aircraft.b(3).vx;
    Payload.vy=aircraft.b(3).vy;
%     disp('Vehicles_1')
%     disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(3).M=aircraft.b(3).M+Payload.M;
    
    % Vehicles   
    L=norm(aircraft.b(4).L);
    l=2;
    nel=(size(aircraft.b(4).M,1)-6)/6;
    Payload=b_constant_p_square(L,l,E,G,rho,nel);
    Payload.o=aircraft.b(4).o;
    Payload.vx=aircraft.b(4).vx;
%     Payload.vy=aircraft.b(4).vy;
%     disp('Vehicles_2')
    disp(sum(sum(Payload.M(1:6:end,1:6:end))))
    aircraft.b(4).M=aircraft.b(4).M+Payload.M;
end

% APU and other systems

E=0;            % Young modulus
G=0;            % Shear modulus
rho=2700/20;    % Density

L=norm(aircraft.b(2).L);
l=1;
nel=(size(aircraft.b(2).M,1)-6)/6;
Payload=b_constant_p_square(L,l,E,G,rho,nel);
Payload.o=aircraft.b(2).o;
Payload.vx=aircraft.b(2).vx;
Payload.vy=aircraft.b(2).vy;
% disp('APU')
% disp(sum(sum(Payload.M(1:6:end,1:6:end))))
aircraft.b(2).M=aircraft.b(2).M+Payload.M;

% Clear unuseful variables
clear human
clear vehicles
clear Payload
clear E
clear G
clear l
clear L
clear nel
clear rho