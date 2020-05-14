% This is a test:
% Transverse vibrations of free-free uniform Euler-Bernoulli beam).
cd ..
cd ..
cd ..

init
clear all, close all, clc
cd Test\'Convergence _studies'\Free_beam\

cd ..
cd ..

%% Data
L=1e3;          % Length of the beam
E=70*1e6;       % Young modulus
G=27*1e6;       % Shear modulus
rho=2700;       % Density alluminium
nel=10;        % Number of elements
l=10;           % Side of the square


%% Exact solution
% load('w_esatta.mat')

cd ..

for i=1:length(nel)
    %% Build model
    tic
    beam=b_constant_p_square(L,l,E,G,rho,nel(i));
    beam.o=[0,0,0]';
    beam.vx=[1,0,0]';
    beam.vy=[0 1 0]';
    beam.oc=true(6,1);
    
    free_beam=m_init;
    node=en_init;
    node.x=[0,0,0]';
    node.K=1e3*sparse(eye(size(node.K)));
    node.M=sparse(zeros(size(node.K)));
    node.C=sparse(zeros(size(node.K)));
    node.c=false(6,1);
    free_beam.en=node;
    
    free_beam=m_add_beam(free_beam,beam);
    
    free_beam=m_compute_matrices(free_beam);
    
    M=free_beam.M;
    K=free_beam.K;
    
    %% Solution
    [w, V, k] = ROM_solver(14, M, K);
    % [V,D,FLAG]=eigs(K,M,size(K,2));
    t(i)=toc;

    %% Save error
    error1=norm(w(1:14)-w_esatta(1:14),1);
    error2=norm(w((1:14)-w_esatta(1:14)),2);
    errorinf=norm(w(1:14)-w_esatta(1:14),'Inf');
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

