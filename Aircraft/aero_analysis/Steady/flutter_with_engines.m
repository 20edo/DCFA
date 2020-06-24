clear all; clc; close all;
cd ..
cd ..
cd generate_model\
inside = 0;
pos_x = [-3.5:-0.25:-7];
pos_z = [-0.5:-0.25:-4];
pos = [];
for i = 1:length(pos_x)
    for j = 1:length(pos_z)
        pos = [pos; pos_x(i), pos_z(j)];
    end
end
pos = [pos; 1,1];
max = size(pos,1)-1;
save pos
cd ..
cd aero_analysis\Steady\
velocity = zeros(max+1,1);
velocity(end) = 1;
save velocity
%% Cycle in order to make:
% discharge_from_inside
% discharge_from_outside
if 1
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
save flutter_engines
%%
clear all
load flutter_engines
pos_x = [-3.5:-0.25:-7];
pos_z = [-0.5:-0.25:-4];
velocity = reshape(velocity(1:end-1),length(pos_z),length(pos_x));
velocity = [pos_z',velocity];
velocity = [velocity; 0, pos_x];
%% Plot
close all
if 1
    figure(1)
    contour(pos_x,pos_z,velocity(1:end-1,2:end),'ShowText','on')
    set(gca,'xdir','reverse')
    axis equal
    contour(pos_x,pos_z,velocity(1:end-1,2:end),'ShowText','on')
    set(gca,'xdir','reverse')
    grid on
    hold on
    p = plot(-5.5,-2,'x','MarkerSize',15,'LineWidth',2,'color',[0.6350, 0.0780, 0.1840])
    xlabel('Horizontal offset $[m]$','fontsize',14,'Interpreter','latex')
    ylabel('Vertical offset $[m]$','fontsize',14,'Interpreter','latex')
    title('h = $10000$ m','fontsize',14,'Interpreter','latex')
    legend(p(1),'Actual position')
    set(gcf, 'Position',  [0, 0, 2000, 2000])
    saveas(figure(1),'flutter_engine','epsc')
end

load handel
sound(y,Fs)