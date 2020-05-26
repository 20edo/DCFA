% just an exercise to compute the divergence and control reversal for a
% straigth wing 
%% Only torsion matters because the beam is straight
clear all 
syms x GJ Le e c CL0 CLa CMac CLb CMb real 
nel = 100; 
 
N = [(1-x)/2 , (1+x)/2]; 
N_x = diff(N)*2/Le; 

% Matrices for the elements 
K_el = int(N_x'*GJ*N_x*Le/2,x,-1,1); 
fa_el = int(N'*(e*c*CL0 + c^2*CMac)*Le/2,x,-1,1);
Ka_el = int(N'*(e*c*CLa)*N*Le/2,x,-1,1); 
fb_el = int(N'*(e*c*CLb + c^2*CMb)*Le/2,x,-1,1); 
Lq_el = int(c*CLa*N*Le/2,x,-1,1); 
Lb = int(c*CLb*Le/2,x,-1,1); 

GJ = 1e7; 
L = 5; 
Le = L/nel; 
e = 0.125; 
c = 1; 
CLa = 2*pi; 
CL0 = 0.1; 
CMac = -pi/2; 
CLb = 2*pi/4; 
CMb = -0.4; 
K_el = eval(K_el); 
fa_el = eval(fa_el); 
Ka_el = eval(Ka_el); 
fb_el = eval(fb_el);
Lq_el = eval(Lq_el); 
Lb = eval(Lb); 
%% Ora devo assemblare le matrici per la beam... 
K = zeros(nel+1,nel+1); 
Ka = zeros(nel+1,nel+1);
fa = zeros(nel+1,1); 
fb = zeros(nel+1,1); 
Lq = zeros(1,nel+1); 
for i = 1:nel 
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ K_el;
    Ka(i:i+1,i:i+1) = Ka(i:i+1,i:i+1)+ Ka_el;
    fa(i:i+1,1) = fa(i:i+1,1) + fa_el;
    fb(i:i+1,1) = fb(i:i+1,1) + fb_el; 
    Lq(1,i:i+1) = Lq(1,i:i+1) + Lq_el;
end

%% Clamp everything 
K = K(2:end,2:end);
Ka = Ka(2:end,2:end);
fa = fa(2:end,1);
fb = fb(2:end,1); 
Lq = Lq(1,2:end);



%% Problem of divergence
[V_div,D_div] = eig(K,Ka); 
lambda_div = sort(diag(D_div)); 
q_div = lambda_div(end-sum(diag(D_div)>1e-6)+1); 
q_ex = pi^2/4*((GJ/L)/(e*L*c*CLa));
disp('q_div = ')
disp(q_div)
disp('q_ex = ')
disp(q_ex)


%% Problem of control reversal
B = [Ka, fb; Lq, Lb]; 
A = zeros(size(B)); 
A(1:size(K,1),1:size(K,2)) = A(1:size(K,1),1:size(K,2)) + K; 

[V_inv,D_inv] = eig(A,B); 
lambda_inv = sort(diag(D_inv)); 
q_inv = lambda_inv(end-sum(diag(D_inv)>1e-6)+1);
disp('q_inv = ')
disp(q_inv)


