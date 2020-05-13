% This is a test:
% Transverse vibrations of clamped-free (cantilever) uniform Euler-Bernoulli beam).
cd ..
cd ..
init
clear all, close all, clc
cd Test\
%% Data
L=1e3;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=round(linspace(4,500,10));        % Number of elements
l=10;           % Side of the square


%% Exact solution
load('w_esatta.mat')
error=zeros(size(nel));

cd ..

for i=1:length(nel)
    %% Build beam
    tic
    beam=b_constant_p_square(L,l,E,G,rho,nel(i));
    %% Constraints
    K=beam.K(7:end,7:end);
    M=beam.M(7:end,7:end);

    %% Solution
    [w, V, k] = ROM_solver(14, M, K);
    % [V,D,FLAG]=eigs(K,M,size(K,2));
    t(i)=toc;

    %% Save error
    error1=norm(w([1,3,5,9,11])-w_esatta([1,3,5,9,11]),1);
    error2=norm(w([1,3,5,9,11])-w_esatta([1,3,5,9,11]));
    errorinf=norm(w([1,3,5,9,11])-w_esatta([1,3,5,9,11]),'Inf');
    e1(i)=error1;
    e2(i)=error2;
    einf(i)=errorinf;
    it(i)=k;
    
end

subplot(2,2,1)
    loglog(nel,einf)
    grid on
    xlabel('Number of elements')
    ylabel('Error')
    title('Norm INF error')
subplot(2,2,2)
    loglog(nel,e2)
    grid on
    xlabel('Number of elements')
    ylabel('Error')
    title('Norm 2 error')
subplot(2,2,3)
    loglog(nel,t)
    grid on
    xlabel('Number of elements')
    ylabel('Time')
    title('Time spent')
subplot(2,2,4)
    plot(nel,it)
    grid on
    xlabel('Number of elements')
    ylabel('ROM_iterations')

