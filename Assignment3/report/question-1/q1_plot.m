close all;

hold on;

bh = csvread('q1data.csv', 1, 0);
plot(bh(1:6, 1), bh(1:6, 2), 'x');

interpol = csvread('q1a.csv', 1, 0);
plot(interpol(:, 1), interpol(:, 2));

xlabel('B (T)');
ylabel('H (A/m)');
legend('Data', 'Interpolation');
grid on;

saveas(gcf, 'q1a.eps', 'epsc');

figure;
hold on;

bh = csvread('q1data.csv', 1, 0);
plot(bh(:, 1), bh(:, 2), 'x');

interpol = csvread('q1b.csv', 1, 0);
plot(interpol(:, 1), interpol(:, 2));

xlabel('B (T)');
ylabel('H (A/m)');
legend('Data', 'Interpolation');
grid on;

saveas(gcf, 'q1b.eps', 'epsc');

figure;
hold on;

bh = csvread('q1data.csv', 1, 0);
plot(bh(:, 1), bh(:, 2), 'x');

interpol = csvread('q1c.csv', 1, 0);
plot(interpol(:, 1), interpol(:, 2));

xlabel('B (T)');
ylabel('H (A/m)');
legend('Data', 'Interpolation');
grid on;

saveas(gcf, 'q1c.eps', 'epsc');
