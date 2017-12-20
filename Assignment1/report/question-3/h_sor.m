clear all; close all; clc;

%% Read data from CSV in the same folder
data = csvread('h_sor.csv', 1, 0);
h = data(:, 1);
iter = data(:, 2);
v = data(:, 3);

%% Generate plot
figure;
plot(h, iter);
xlabel('1/h');
ylabel('Iterations');
grid on;

saveas(gcf, 'h_iter_sor.eps', 'epsc');

%% Generate plot
figure;
plot(h, v);
xlabel('1/h');
ylabel('V (V)');
grid on;

saveas(gcf, 'h_v_sor.eps', 'epsc');

