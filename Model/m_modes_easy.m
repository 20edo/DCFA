function m_modes_easy(model,N,alpha)

if (~exist('alpha', 'var'))
    % "dof" parameter does not exist
    alpha = 0.3;
end

[U, lambda] = eigs(model.K+alpha*model.M,model.M,N,'smallestabs');
w = sqrt(diag(lambda)-alpha);
for i = 1:length(model.b)
    u_beam = transpose(model.b(i).A)*U*100;
    for k = 1:model.b(i).nel+1
        model.b(i).in(k).d = u_beam(1+6*(k-1):6*(k),:);
    end
end
clear i k
for j = N
    figure
    for i=1:length(model.b)
        n = length(model.b(i).in);
        for k = 1:n
            temp1 = model.b(i).o + model.b(i).vx*model.b(i).in(k).x;
            plot3(temp1(1),temp1(2),temp1(3),'o','Color',[0, 0.4470, 0.7410]);
            hold on
            temp = model.b(i).o + model.b(i).vx*model.b(i).in(k).x + model.b(i).in(k).d(1:3,j);
            plot3(temp(1),temp(2),temp(3),'o','Color',[0.8500, 0.3250, 0.0980]);
        end
    end
    hold off
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
end