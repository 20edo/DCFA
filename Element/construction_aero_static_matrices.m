clear all; close all; clc; 
syms Lambda x GJ EJ Le e c CL0 CLa CMac CLb CMb real
N_th = [(1-x)/2 , (1+x)/2]; 
N_w = [1/4*(2-3*x+x^3),-1/4*(1-x-x^2+x^3)*Le/2,1/4*(2+3*x-x^3),-1/4*(-1-x+x^2+x^3)*Le/2]; 

K_thth = int(diff(N_th)'*2/Le*GJ*diff(N_th)*2/Le,x,-1,1)*Le/2; 
K_ww = int(diff(diff(N_w'))*4/Le^2*EJ*diff(diff(N_w))*4/Le^2,x,-1,1)*Le/2; 

%% Right wing
fa_th = int(N_th'*cos(Lambda)^2*(c^2*CMac + e*c*CL0),x,-1,1)*Le/2; 
fa_w = -int(N_w'*cos(Lambda)*c*CL0,x,-1,1)*Le/2 + ...
    int(diff(N_w')*2/Le*cos(Lambda)*sin(Lambda)*(c^2*CMac + e*c*CL0),x,-1,1)*Le/2; 
Ka_thth = int(N_th'*cos(Lambda)^3*e*c*CLa*N_th,x,-1,1)*Le/2; 
Ka_thw = int(N_th'*cos(Lambda)^2*sin(Lambda)*e*c*CLa*diff(N_w)*2/Le,x,-1,1)*Le/2; 
Ka_wth = -int(N_w'*cos(Lambda)^2*c*CLa*N_th,x,-1,1)*Le/2 + ... 
    int(diff(N_w')*2/Le*cos(Lambda)^2*sin(Lambda)*c*CLa*N_th,x,-1,1)*Le/2; 
Ka_ww = -int(N_w'*cos(Lambda)*sin(Lambda)*c*CLa*diff(N_w)*2/Le,x,-1,1)*Le/2 + ...
    int(diff(N_w')*2/Le*cos(Lambda)*sin(Lambda)^2*c*CLa*diff(N_w)*2/Le,x,-1,1)*Le/2; 