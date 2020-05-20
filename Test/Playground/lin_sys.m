function sist = lin_sys(t,y,M,C,K,f)
n2 = size(M,1); 
sist = zeros(2*n2,1);
sist(1:n2) = y(n2+1:end); 
M = diag(M);
C = diag(C); 
K = diag(K); 
sist(n2+1:end) = 1./M(1:end).*(-C(1:end).*y(n2+1:end)-K(1:end).*y(1:n2)+f(t)); 
end