function m_plot_easy(model) 

for j = 1:20:size(model.b(1).in(1).d,2)
    figure(1)
for i=1:length(model.b)
    n = length(model.b(i).in); 
    for k = 1:n
        
       temp1 = model.b(i).o + model.b(i).vx*model.b(i).in(k).x;
       plot3(temp1(1),temp1(2),temp1(3),'o','Color',[0, 0.4470, 0.7410]); 
        hold on 
        xlim([-200,200]);
        ylim([-220,220]);
        zlim([-300,500]);
        view(50,20);
       temp = model.b(i).o + model.b(i).vx*model.b(i).in(k).x + model.b(i).in(k).d(1:3,j);
       plot3(temp(1),temp(2),temp(3),'o','Color',[0.8500, 0.3250, 0.0980]);   
    end
end
hold off
disp(j)
end

