close all;

data = csvread('q3b.csv', 1, 0);
iteration = data(:, 1);
v1 = data(:, 2);
v2 = data(:, 3);
f1 = data(:, 4);
f2 = data(:, 5);
error = data(:, 6);

hold on;

plot(iteration, v1-v2, iteration, v2);
xlabel('Iteration');
ylabel('Voltage (V)');
legend('V_A', 'V_B');
grid on;

saveas(gcf, 'voltages.eps', 'epsc');

figure;
hold on;

plot(iteration, f1, iteration, f2);
xlabel('Iteration');
ylabel('f(v)');
legend('f_1', 'f_2');
grid on;

saveas(gcf, 'residual.eps', 'epsc');

figure;
hold on;

plot(iteration, error);
xlabel('Iteration');
ylabel('Error');
grid on;

saveas(gcf, 'error.eps', 'epsc');

figure;
hold on;

plot(iteration(2:end), error(2:end)./error(1:end-1).^2);
xlabel('Iteration');
ylabel('\epsilon_k / \epsilon_{k-1}^2');
grid on;

saveas(gcf, 'error_ratio.eps', 'epsc');
