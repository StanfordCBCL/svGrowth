% Given values
t0 = 2;
p0 = 93;
t1 = 10;
deltap = 50;
P_h = 100
delta_P= 1+0.5;

p_slope = (P_h*delta_P -P_h) / (t1 - 2); 
% Define the linear gradual loading equation
p = @(t) P_h + p_slope * (t - 2);

% Time vector for plotting
t = linspace(t0 - 2, t1 + 2, 100);

% Plot
plot(t, p(t), 'b-', 'LineWidth', 2);
hold on;
plot([t0 t1], [p0 p0*delta_P], 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('t');
ylabel('p(t)');
title('Line passing through (t0, p0) and (t1, p0 + deltap)');
grid on;
