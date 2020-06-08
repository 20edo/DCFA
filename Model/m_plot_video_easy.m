function FR=m_plot_video_easy(model,t,U)
%% Add the zeros for the constrained displacements
for j=1:length(model.en)
    for k=1:6
        if model.en(j).c(k) 
            index = 6*(j-1)+k;
            U = [U(1:index-1,:);zeros(1,size(U,2));U(index:end,:)]; 
        end
    end            
end
%% Add the eigenshapes into the displacements
for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*U; 
    for k = 1:model.b(i).nel+1
       model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end    
end
clear i k
h = length(t);
% FR(k) = struct('cdata',[],'colormap',[]);
for j = 1:h
    temp = [];
    temp1 = [];
    for i=1:length(model.b)
        n = length(model.b(i).in);
        for k = 1:n
            temp1 =[temp1 model.b(i).o + model.b(i).vx*model.b(i).in(k).x];
            temp = [temp model.b(i).o + model.b(i).vx*model.b(i).in(k).x + model.b(i).in(k).d(1:3,j)];
        end
    end
    plot3(temp1(1,:),temp1(2,:),temp1(3,:),'o','Color',[0, 0.4470, 0.7410]);
    hold on
    plot3(temp(1,:),temp(2,:),temp(3,:),'o','Color',[0.8500, 0.3250, 0.0980]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title({['Deformed model at time:',num2str(t(j))]});
    FR(j) = getframe;
    hold off
end
delta_t=t(2)-t(1);
movie(FR,1,1/delta_t)
