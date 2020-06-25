clear all; clc; close all;
cd ..
cd ..
cd generate_model\
%% Problem
% 1 -> discharge from inside
% 2 -> discharge from outside
% 3 -> discharge all toghether
problem = 3;
%%%% !!! CHANGE PROBLEM ALSO DOWN BEFORE RUNNING !!! %%%%
switch problem
    case 1
        carb = [    1.0000    1.0000    1.0000
            0.9000    1.0000    1.0000
            0.8000    1.0000    1.0000
            0.7000    1.0000    1.0000
            0.6000    1.0000    1.0000
            0.5000    1.0000    1.0000
            0.4000    1.0000    1.0000
            0.3000    1.0000    1.0000
            0.2000    1.0000    1.0000
            0.1000    1.0000    1.0000
            0    1.0000    1.0000
            0    0.9000    1.0000
            0    0.8000    1.0000
            0    0.7000    1.0000
            0    0.6000    1.0000
            0    0.5000    1.0000
            0    0.4000    1.0000
            0    0.3000    1.0000
            0    0.2000    1.0000
            0    0.1000    1.0000
            0         0    1.0000
            0         0    0.9000
            0         0    0.8000
            0         0    0.7000
            0         0    0.6000
            0         0    0.5000
            0         0    0.4000
            0         0    0.3000
            0         0    0.2000
            0         0    0.1000
            0         0         0
            1         1         1];
    case 2
        carb = [    1.0000    1.0000    1.0000
            0.9000    1.0000    1.0000
            0.8000    1.0000    1.0000
            0.7000    1.0000    1.0000
            0.6000    1.0000    1.0000
            0.5000    1.0000    1.0000
            0.4000    1.0000    1.0000
            0.3000    1.0000    1.0000
            0.2000    1.0000    1.0000
            0.1000    1.0000    1.0000
            0    1.0000    1.0000
            0    0.9000    1.0000
            0    0.8000    1.0000
            0    0.7000    1.0000
            0    0.6000    1.0000
            0    0.5000    1.0000
            0    0.4000    1.0000
            0    0.3000    1.0000
            0    0.2000    1.0000
            0    0.1000    1.0000
            0         0    1.0000
            0         0    0.9000
            0         0    0.8000
            0         0    0.7000
            0         0    0.6000
            0         0    0.5000
            0         0    0.4000
            0         0    0.3000
            0         0    0.2000
            0         0    0.1000
            0         0         0
            1         1         1];
        carb = flip(carb,2);
    case 3
        temp = linspace(1,0,31);
        carb(:,1) = temp';
        carb(:,2) = temp';
        carb(:,3) = temp';
        carb = [carb; 1,1,1]; 
end

max = size(carb,1)-1;
save carb
cd ..
cd aero_analysis\Steady\
velocity = zeros(max+1,1);
velocity(end) = 1;
save velocity
%% Cycle in order to make:
% discharge_from_inside
% discharge_from_outside
if 0
    contator = 1;
    while contator<=length(velocity)-1
        flutter_study_as_linear
        flutter_velocity = v(i);
        load velocity
        contator = velocity(end)
        velocity(contator) = flutter_velocity;
        contator = contator + 1;
        velocity(end) = contator;
        clearvars -except velocity contator
        save velocity
    end
    problem = 3; 
    switch problem
        case 1
            save discharge_from_inside.mat
        case 2
            save discharge_from_outside.mat
        case 3
            save discharge_all.mat
    end
end

%% Display results
close all 
fuel_max = [2.2687e+04; 8.7519e+03; 1.0802e+04];
carb = [    1.0000    1.0000    1.0000
    0.9000    1.0000    1.0000
    0.8000    1.0000    1.0000
    0.7000    1.0000    1.0000
    0.6000    1.0000    1.0000
    0.5000    1.0000    1.0000
    0.4000    1.0000    1.0000
    0.3000    1.0000    1.0000
    0.2000    1.0000    1.0000
    0.1000    1.0000    1.0000
    0    1.0000    1.0000
    0    0.9000    1.0000
    0    0.8000    1.0000
    0    0.7000    1.0000
    0    0.6000    1.0000
    0    0.5000    1.0000
    0    0.4000    1.0000
    0    0.3000    1.0000
    0    0.2000    1.0000
    0    0.1000    1.0000
    0         0    1.0000
    0         0    0.9000
    0         0    0.8000
    0         0    0.7000
    0         0    0.6000
    0         0    0.5000
    0         0    0.4000
    0         0    0.3000
    0         0    0.2000
    0         0    0.1000
    0         0         0
    1         1         1];
load discharge_from_inside.mat
for i=1:length(velocity)-1
    perc(i) = carb(i,:)*fuel_max/(sum(fuel_max))*100;
end
perc_in = 100-perc;
velocity_in = velocity(1:end-1);
% figure
% plot(perc_in,velocity_in)


carb = [    1.0000    1.0000    1.0000
    0.9000    1.0000    1.0000
    0.8000    1.0000    1.0000
    0.7000    1.0000    1.0000
    0.6000    1.0000    1.0000
    0.5000    1.0000    1.0000
    0.4000    1.0000    1.0000
    0.3000    1.0000    1.0000
    0.2000    1.0000    1.0000
    0.1000    1.0000    1.0000
    0    1.0000    1.0000
    0    0.9000    1.0000
    0    0.8000    1.0000
    0    0.7000    1.0000
    0    0.6000    1.0000
    0    0.5000    1.0000
    0    0.4000    1.0000
    0    0.3000    1.0000
    0    0.2000    1.0000
    0    0.1000    1.0000
    0         0    1.0000
    0         0    0.9000
    0         0    0.8000
    0         0    0.7000
    0         0    0.6000
    0         0    0.5000
    0         0    0.4000
    0         0    0.3000
    0         0    0.2000
    0         0    0.1000
    0         0         0
    1         1         1];

carb = flip(carb,2);

load discharge_from_outside.mat
for i=1:length(velocity)-1
    perc(i) = carb(i,:)*fuel_max/(sum(fuel_max))*100;
end
perc_out = 100-perc;
velocity_out = velocity(1:end-1);
% figure
% plot(perc_out,velocity_out)
clear carb
temp = linspace(1,0,31);
        carb(:,1) = temp';
        carb(:,2) = temp';
        carb(:,3) = temp';
        carb = [carb; 1,1,1]; 
        load discharge_all.mat
for i=1:length(velocity)-1
    perc(i) = carb(i,:)*fuel_max/(sum(fuel_max))*100;
end
perc_all = 100-perc;
velocity_all = velocity(1:end-1);

figure(1)
plot(perc_in,velocity_in,perc_out,velocity_out,perc_all,velocity_all,'LineWidth',1.5)
hold on 
plot(perc_out(8),velocity_out(8),'.','MarkerSize',20,'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
plot(perc_out(11),velocity_out(11),'.','MarkerSize',20,'LineWidth',2,'Color',[0.3010, 0.7450, 0.9330])
ylabel('Flutter velocity \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
xlabel('Fuel discharged \quad $[\%]$','fontsize',14,'interpreter','latex')
title('h = $10000$ m','fontsize',14,'interpreter','latex');
grid on
legend('Discharge T1 $\rightarrow$ T2 $\rightarrow$ T3','Discharge T3 $\rightarrow$ T2 $\rightarrow$ T1','Discharge all tanks together','Condition A','Condition B','fontsize',12,'interpreter','latex')
set(gcf, 'Position',  [0, 0, 800, 500])
saveas(figure(1),'flutter_with_fuel','epsc')