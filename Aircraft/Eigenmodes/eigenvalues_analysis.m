% This script plots the first 200 eigenvalues of the aircraft in
% semi-logaritmic scale

%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%         
%

%% Generate the model
cd ..
cd generate_model\
generate_model
cd Eigenmodes\

%% Calculate frequencies
number=200;     % Number of frequencies required
% Solve
alpha=1;        % Shift
[V,D,flag] = eigs(aircraft.K+alpha*aircraft.M,aircraft.M,number,'smallestabs');
w=real(diag(D-alpha).^0.5);

fig=figure;
set(gcf, 'Position',  [100, 100, 5000, 4000])
semilogy(1:length(w)-6,w(7:end),'xr')
title('Aircraft')
xlabel('Number of eigenvalue')
ylabel('Frequency')
grid on
saveas(fig,'Eigenvalues_aircraft','svg')