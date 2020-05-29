% This script plots the eigenvalues in semilogaritmic scale

clear all, close all, clc

load w_esatta.mat

semilogy(1:length(w_esatta),w_esatta,'xr')
xlabel('Number of eigenvalue')
ylabel('Frequency')

grid on
