% Blasius Solver without a Function

tic;

del_y = 0.0005;
x = 0.5;
Re = 10000;

del = 5 / sqrt(Re);
y = 0:del_y:2 * del;
eta_0 = 0;
eta_max = max(y) * sqrt(Re / x);
fgh = zeros(3, 1);
h = [0.1, 0.3];
d_eta = del_y * sqrt(Re / x);
N = (eta_max - eta_0) / d_eta;
err = 1;

% Initializing solution
eta = eta_0 : d_eta : eta_max;

while err > 1e-15
    for j = 1:2
        F = zeros(N + 1, 3);
        F(1, :) = [0 , 0, h(j)];

        % Trail using RK-4 Method
        for i = 1:N
            k = d_eta * [F(i, 2), F(i, 3), -(F(i, 1) * F(i, 3)) / 2];
            F(i + 1, :) = F(i, :) + (1 / 6) * (k + 2 * k + 2 * k + k);
        end
        g_inf(j) = F(end, 2);
    end

    [err, index] = max(abs(1 - g_inf));
    if diff(g_inf) ~= 0
        h(index) = h(1) + (diff(h) / diff(g_inf)) * (1 - g_inf(1));
    end
end

% Plotting
figure;
plot(y, F(:, 2), 'LineWidth', 1, 'Color', 'b', 'Marker', 'o', 'MarkerSize', 2);
hold on;
plot(y, F(:,3), 'LineWidth', 1, 'Color', 'r', 'Marker', 's', 'MarkerSize', 2);
xlabel('y');
ylabel('Velocity');
title('Blasius Solution');
legend('u', 'v');
grid on;

toc;


