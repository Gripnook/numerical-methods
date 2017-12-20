clear all; close all; clc;

%% Read data from CSV in the same folder
data = csvread('w_sor.csv', 1, 0);
w = data(1:end-1, 1);
iter = data(1:end-1, 2);

%% Generate plot
plot(w, iter);
xlabel('w');
ylabel('Iterations');
grid on;

saveas(gcf, 'w_sor.eps', 'epsc');

