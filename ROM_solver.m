function [w, V, k] = ROM_solver(N, M, K, alpha)

% This function gives the eigenshapes relative to the firsts N eigenvalues    
% given:
% N         is the number of eigenshpaes we want to study;
% M         is the complete Mass matrix;
% K         is the complete Stiffness matrix;
% nel       is the number of element considered;
% L         is the length of the beam;
% alpha     is the shift parameter
% dof       are the number of degree of freedom we wanted to analyse. (if not specified
% it's considered equal to 6: u v w th phi psi);
% plot_cond is a bool condition to decide whether plotting graphs is needed
%
% (to implement: given the structure beam, add lumped mass as engine,...)
%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
%% Check default


%  if (~exist('dof', 'var'))
%      % "dof" parameter does not exist
%       dof = 6;
%  end
%  if (~exist('plot_cond', 'var'))
%      % "plot" parameter does not exist
%       plot_cond = 0;
%  end
 
 if (~exist('alpha', 'var'))
     % "dof" parameter does not exist
      alpha = 0;
 end
 
N_req = N; 
N = max([2*N,N+20]); 
 
%% subspace iter
n = size(K, 1);
V0 = randn(n, N); 
V0 = V0*diag(1./(diag(V0'*V0).^0.5));

Vk = V0;
lambda_k = zeros(N, 1);
Lambda_k = [];

tolerance = 1e-6; 

Lower = chol(K + alpha*M, 'lower');
for k = 1:1000
    % Uk = K \ M * Uk
    % Uk = (L * L') \ M * Uk
    Vk = (transpose(Lower)\(Lower\(M * Vk)));
    Vk_old = Vk; 
    % reduced matrices
    kk = Vk'*K*Vk;
    mm = Vk'*M*Vk;
    
    lk_prev = lambda_k;
    
    [Q, lambda_k] = eig(kk, mm);
    lambda_k = diag(lambda_k);
    [dmy, II] = sort(lambda_k);
    lambda_k = lambda_k(II);
    Q = Q(:, II);
    Lambda_k = [Lambda_k, lambda_k];
    
    Vk = Vk*Q;
    Vk = Vk*diag(1./(diag(Vk'*Vk).^0.5));
    
    if (k == 1)
        % skip first iteration
        continue;
        
    elseif (norm(lambda_k - lk_prev) < tolerance && ...
            max(max(abs(kk-diag(diag(kk)))))/max(diag(kk))< tolerance*1e-5)
        disp(sprintf('ROM solver converged in %d iterations', k));
        break;
    end
end

w = sqrt(lambda_k); 
V = Vk; %eigenvector matrix

if max(abs(imag(w(7:end))))/max(abs(real(w(7:end))))<1e-5
w = real(w(1:N_req)); 
V = V(:,1:N_req); 
else 
    error('Check the eingenvalues')
end
if (k==1000)
     disp(sprintf('ROM solver did not converge'));
end

% %% plot of the dof
% % CHECK: to implement for a 3-D structure
% %this is an example for a beam
% 
% if plot_cond == 1
%     
%     x = [1:nel]/nel*L;
% 
%     %u v w th phi psi
%     for i = 1:N
%         figure;
%         subplot(2, 3, 1);
%         plot(x, V(1:dof:end, i));
%         xlabel('x');
%         ylabel('u');
%         title(sprintf('mode %d: \\omega=%e', i, w(i)));
%     
%         subplot(2, 3, 2);
%         plot(x, V(2:dof:end, i));
%         xlabel('x');
%         ylabel('v');
%     
%         subplot(2, 3, 3);
%         plot(x, V(3:dof:end, i));
%         xlabel('x');
%         ylabel('w');
%     
%         subplot(2, 3, 4);
%         plot(x, V(4:dof:end, i));
%         xlabel('x');
%         ylabel('\theta');
%     
%         subplot(2, 3, 5);
%         plot(x, V(5:dof:end, i));
%         xlabel('x');
%         ylabel('\phi');
%     
%         subplot(2, 3, 6);
%         plot(x, V(6:dof:end, i));
%         xlabel('x');
%         ylabel('\psi');
%     
%     end
% end



