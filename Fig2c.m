function I1_FFL_Simulation_Hill
% parameters
k1 = 1;
k2 = 1;
k3 = 1;
k5 = 1;
n = 1;

% initial conditions
X0 = 1; % concentration of X
Y0 = 0; % initial concentration of Y
Z0 = 0; % initial concentration of Z

tspan = [0 16]; % time span for the simulation
init = [Y0; Z0]; % initial condition vector

% Define ODEs and Solve for each value of k4
k4_values = [1.5,0.5,0.01];
colors = {[1 0 0], [1 0.5 0.5], [1 0.65 0.65], [1 0.8 0.8]}; % cell array of color specifications


figure('Position', [100, 100, 500,375]);
hold on
for i = 1:length(k4_values)
    k4 = k4_values(i);
    odes = @(t, y) [k1*X0 - k2*y(1); k3*X0/(1 + (y(1)/k4)^n) - k5*y(2)];
    [t, y] = ode45(odes, tspan, init);
    % Normalize the output concentration (Z)
    y(:,2) = (y(:,2) - min(y(:,2))) / (max(y(:,2)) - min(y(:,2)));
    % Plot results for each value of k4
    plot(t, y(:, 2), 'Color', colors{i}, 'LineWidth', 2)
    k4_labels{i} = sprintf('IFFL (k_n = %.2f)', k4);
end
% Add simple regulation model
odes_simple = @(t, Z) k3*X0 - k5*Z;
[t_simple, Z_simple] = ode45(odes_simple, tspan, Z0);
Z_simple = (Z_simple - min(Z_simple)) / (max(Z_simple) - min(Z_simple)); % Normalize
plot(t_simple, Z_simple, 'k--', 'LineWidth', 2) % Plot with dashed black line

ylim([0 1.1]) % adjust the y-axis limits
legend([k4_labels, 'Direct Regulation'], 'Location', 'best')
xlabel('Cell Generations')
xlim([-1 5])
ylabel('Normalized Z')
%title('Incoherent Type I Feed-Forward Loop vs. Direct Regulation')
%title('IFFL vs. Direct Regulation')
%box on
%set(gca, 'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
end
