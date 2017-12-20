clear all; close all; clc;

%% Read data from CSV in the same folder
data = csvread('banded.csv', 1, 0);
N = data(:, 1);
R = data(:, 2);
ff = fittype('a * log(x + c) + b');
f = fit(N, R, ff);

%% Generate plot
plot(f, N, R);
xlabel('N');
ylabel('R (k\Omega)');
legend('data points', strcat(string(f.a), ' log(N+', string(f.c), ')+', string(f.b)), 'Location', 'northwest');
grid on;

saveas(gcf, 'resistance.eps', 'epsc');

