close all;

data = csvread('q4a.csv', 1, 0);
segments = data(:, 1);
error = data(:, 2);

disp(mean(diff(log(error)) ./ diff(log(segments))));

hold on;

plot(log(segments), log(error));
xlabel('ln(N)');
ylabel('ln(E)');
grid on;

saveas(gcf, 'q4a.eps', 'epsc');

figure;

data = csvread('q4b.csv', 1, 0);
segments = data(:, 1);
error = data(:, 2);

disp(mean(diff(log(error)) ./ diff(log(segments))));

hold on;

plot(log(segments), log(error));
xlabel('ln(N)');
ylabel('ln(E)');
grid on;

saveas(gcf, 'q4b.eps', 'epsc');

figure;

data = csvread('q4c.csv', 1, 0);
r = data(:, 1);
error = data(:, 3);

hold on;

plot(r, error);
xlabel('r');
ylabel('E');
grid on;

saveas(gcf, 'q4c.eps', 'epsc');
