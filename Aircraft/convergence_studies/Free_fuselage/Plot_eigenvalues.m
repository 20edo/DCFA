% This script plots the eigenvalues in semilogaritmic scale

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

clear all, close all, clc

load w_esatta.mat

fig=figure;
semilogy(1:length(w_esatta)+6,[0 0 0 0 0 0 w_esatta(1:end)'],'xr')
xlabel('Number of eigenvalue')
ylabel('Frequency')
grid on
saveas(fig,'Eigenvalues_free_fuselage','epsc')