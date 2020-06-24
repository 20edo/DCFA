clear all; clc; close all;
cd ..
cd ..
cd generate_model\
inside = 0;
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
save carb
if ~inside
    carb = flip(carb,2);
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
        flutter_study2
        flutter_velocity = v(i);
        load velocity
        contator = velocity(end)
        velocity(contator) = flutter_velocity;
        contator = contator + 1;
        velocity(end) = contator;
        clearvars -except velocity contator
        save velocity
    end
end

%% Display results
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

plot(perc_in,velocity_in,perc_out,velocity_out,'LineWidth',1.5)
ylabel('Flutter velocity \quad $[\frac{m}{s}]$','fontsize',14,'interpreter','latex')
xlabel('Fuel discharged \quad $[\%]$','fontsize',14,'interpreter','latex')
title('h = $10000$ m','fontsize',14,'interpreter','latex');
grid on
legend('Discharged from inside','Discharged from outside','interpreter','latex')
set(gcf, 'Position',  [0, 0, 700, 400])