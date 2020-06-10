function res=funz(m,v,V,q,scaling,X,e)

%% Compute loads
chord=7.72;
k=chord/2*imag(e)/norm(v);
m=m_add_unsteady_loads(m,[1 0 0]',k);

%% Reduce matrices
alpha = 0;
gamma = 0;
Cs = alpha*m.M + gamma*m.K;

M = V'*m.M*V/scaling;
K = V'*m.K*V/scaling;
Cs = V'*Cs*V/scaling;
Ham=V'*m.Ham*V/scaling;
%% Calculate residual
res=(e^2*M+e*Cs+K-q*Ham)*X;

end