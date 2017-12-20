clear all; close all; clc;

%% Read data from CSV in the same folder
data = csvread('cg_results.csv', 1, 0);
iter = data(:, 1);
twoNorm = data(:, 2);
infNorm = data(:, 3);

%% Generate plot
plot(iter, twoNorm, iter, infNorm);
xlabel('Iterations');
ylabel('Residual');
legend('2-Norm', 'Infinity Norm');
grid on;

saveas(gcf, 'cg_results.eps', 'epsc');
