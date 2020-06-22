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
f=w/2/pi;
fig=figure;
set(gcf, 'Position',  [0, 0, 600, 400])
semilogy(1:length(f)-6,f(7:end),'x')
% title('Aircraft')
xlabel('Number of eigenvalue','Interpreter','latex')
ylabel('Frequency \quad [Hz]','Interpreter','latex')
grid on
xlim([0 100])
saveas(fig,'Eigenvalues_aircraft','epsc')