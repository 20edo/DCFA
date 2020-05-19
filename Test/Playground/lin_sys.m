function sist = lin_sys(t,y,M,C,K,f)
n = size(M,1); 
sist = zeros(2*n,1);
sist(1:n) = y(n+1:end); 
M = diag(M);
C = diag(C); 
K = diag(K); 
sist(n+1:end) = 1./M(1:end).*(-C(1:end).*y(n+1:end)-K(1:end).*y(1:n)+f(t)); 