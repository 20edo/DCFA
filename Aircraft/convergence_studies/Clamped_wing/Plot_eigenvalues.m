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
semilogy(1:length(w_esatta),w_esatta,'xr')
title('Clamped wing')
xlabel('Number of eigenvalue')
ylabel('Frequency')
grid on
saveas(fig,'Eigenvalues_clamped_wing','svg')