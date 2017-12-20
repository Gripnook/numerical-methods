clear all; close all; clc;

%% Read data from CSV in the same folder
data = csvread('standard.csv', 1, 0);
N = data(:, 1);
time = data(:, 3);
ff = fittype('a * x ^ b');
f = fit(N, time, ff);

%% Generate plot
plot(f, N, time);
xlabel('N');
ylabel('Time (us)');
legend('data points', strcat(string(f.a), ' N^{', string(f.b), '}'), 'Location', 'northwest');
grid on;

saveas(gcf, 'standard_time.eps', 'epsc');

