function I1_FFL_Fitting
% Load and preprocess experimental data
filename = 'SupFig3_Fitting.xlsx';
sheet = 1;
range = 'A3:G121';
exp_data = xlsread(filename, sheet, range);

% Convert time from minutes to cell generations
t_exp = exp_data(:, 1) / 40;
data_direct = exp_data(:, 2:4);
data_I1FFL = exp_data(:, 5:7);

% Filter data to the range of 1 to 15 cell generations
indices = (t_exp >= 1) & (t_exp <= 15);
t_exp = t_exp(indices);
data_direct = data_direct(indices, :);
data_I1FFL = data_I1FFL(indices, :);

% Compute the average and standard deviation of the Direct Regulation data
avg_direct = mean(data_direct, 2);
std_direct = std(data_direct, 0, 2);

% Normalizing the Direct Regulation data
min_direct = min(avg_direct);
max_direct = max(avg_direct);
avg_direct_normalized = (avg_direct - min_direct) / (max_direct - min_direct);
std_direct_normalized = std_direct / (max_direct - min_direct);

% Compute the average and standard deviation of the I1-FFL data
avg_I1FFL = mean(data_I1FFL, 2);
std_I1FFL = std(data_I1FFL, 0, 2);

% Normalizing the I1-FFL data
min_I1FFL = min(avg_I1FFL);
max_I1FFL = max(avg_I1FFL);
avg_I1FFL_normalized = (avg_I1FFL - min_I1FFL) / (max_I1FFL - min_I1FFL);
std_I1FFL_normalized = std_I1FFL / (max_I1FFL - min_I1FFL);

% Fit curves for Direct Regulation and I1-FFL
params_direct = fit_parameters_direct(t_exp, avg_direct_normalized);
params_I1FFL = fit_parameters_I1FFL(t_exp, avg_I1FFL_normalized);

% Simulate the fitted curves
t_fit = linspace(min(t_exp), max(t_exp), 100);
y_direct_fit = ode_simulation_direct(params_direct, t_fit);
y_I1FFL_fit = ode_simulation_I1FFL(params_I1FFL, t_fit);

% Print out the fitting parameters
fprintf('Fitting parameters for Direct Regulation: k3 = %.3f, n = %.3f, k5 = %.3f\n', params_direct(1), params_direct(2), params_direct(3));
fprintf('Fitting parameters for I1-FFL: k1 = %.3f, k2 = %.3f, k3 = %.3f, k4 = %.3f, k5 = %.3f, n = %.3f\n', params_I1FFL(1), params_I1FFL(2), params_I1FFL(3), params_I1FFL(4), params_I1FFL(5), params_I1FFL(6));

% Plotting
% Plot the results
figure('Position', [100, 100, 500, 375]);
hold on

% Plot average line for Direct Regulation
plot(t_exp, avg_direct_normalized, 'k', 'LineWidth', 2)

% Plot average line for I1-FFL
plot(t_exp, avg_I1FFL_normalized, 'r', 'LineWidth', 2)

% Plot fitted curves for Direct Regulation and I1-FFL
plot(t_fit, y_direct_fit, 'k--', 'LineWidth', 2)
plot(t_fit, y_I1FFL_fit, 'r--', 'LineWidth', 2)

% Plot 50% of average line for I1-FFL
yline(0.455, 'LineWidth', 2);
hold off
bgColor = [0.68 0.85 0.9 0.5]; % Light blue color
set(gca, 'Color', bgColor);

% Add labels and title
xlabel('Cell Generations')
ylabel('Z/Z_{st}')
legend('Direct Regulation', 'IFFL', 'Direct Regulation Fitting' ,'IFFL Fitting', '50% Steady-State of IFFL','Location', 'southeast')
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold')
box on

% Adjust the y-axis limits if needed
ylim([0, 1.1])

y_target = 0.455;  % y-value we are interested in
% Find the nearest y-values above and below the target value
idx_below_1 = find(y_I1FFL_fit <= y_target, 1, 'last');
idx_above_1 = find(y_I1FFL_fit >= y_target, 1, 'first');
idx_below_2 = find(y_direct_fit <= y_target, 1, 'last');
idx_above_2 = find(y_direct_fit >= y_target, 1, 'first');

% Interpolate x at the target y-value between these two points
x_target_I1FFl = interp1([y_I1FFL_fit(idx_below_1), y_I1FFL_fit(idx_above_1)], [t_fit(idx_below_1), t_fit(idx_above_1)], y_target);
fprintf('For y_I1FFL_fit = 0.455, the corresponding x (time) value is: %.3f cell generations\n', x_target_I1FFl);
x_target_Dir = interp1([y_direct_fit(idx_below_2), y_direct_fit(idx_above_2)], [t_fit(idx_below_2), t_fit(idx_above_2)], y_target);
fprintf('For y_direct_fit = 0.455, the corresponding x (time) value is: %.3f cell generations\n', x_target_Dir);
end

function params_direct = fit_parameters_direct(t_exp, avg_direct)
% Initial guess for parameters
k3_guess = 1;
n_guess = 2;
k5_guess = 1;

% Set lower and upper bounds for parameters
lb = [0, 1, 0];
ub = [Inf, 5, Inf];

% Objective function for parameter estimation
objective_func = @(params) compute_error_direct(params, t_exp, avg_direct);

% Fit parameters using lsqcurvefit
params_direct = lsqcurvefit(@(params,t) ode_simulation_direct(params, t), [k3_guess, n_guess, k5_guess], t_exp, avg_direct, lb, ub);

    function error = compute_error_direct(params, t, data)
        k3 = params(1);
        n = params(2);
        k5 = params(3);
        X0 = 2e-6;  % assuming a constant X0
        y_sim = ode_simulation_direct([k3, n, k5, X0], t);
        error = y_sim - data;
    end
end

function params_I1FFL = fit_parameters_I1FFL(t_exp, avg_I1FFL)
% Initial guess for parameters
k1_guess = 1;
k2_guess = 1;
k3_guess =10;
k4_guess =300;
k5_guess =10;
n_guess = 2;

% Set lower and upper bounds for parameters
lb = [0, 0, 0, 0, 0, 1];
ub = [Inf, Inf, Inf, Inf, Inf, 5];

% Objective function for parameter estimation
objective_func = @(params) compute_error_I1FFL(params, t_exp, avg_I1FFL);

% Fit parameters using lsqcurvefit
params_I1FFL = lsqcurvefit(@(params,t) ode_simulation_I1FFL(params, t), [k1_guess, k2_guess, k3_guess, k4_guess, k5_guess, n_guess], t_exp, avg_I1FFL, lb, ub);

    function error = compute_error_I1FFL(params, t, data)
        k1 = params(1);
        k2 = params(2);
        k3 = params(3);
        k4 = params(4);
        k5 = params(5);
        n = params(6);
        X0 = 2e-6;
        y_sim = ode_simulation_I1FFL([k1, k2, k3, k4, k5, n, X0], t);
        error = y_sim - data;
    end
end

function y_sim_direct = ode_simulation_direct(params, t_sim)
k3 = params(1);
n = params(2);
k5 = params(3);
X0 = 2e-6;

% Define ODEs
odefun = @(t,y) k3*X0 - k5*y;

% Solve ODEs
[~, y_sim_direct] = ode45(odefun, t_sim, 0);
end

function y_sim_I1FFL = ode_simulation_I1FFL(params, t_sim)
k1 = params(1);
k2 = params(2);
k3 = params(3);
k4 = params(4);
k5 = params(5);
n = params(6);
X0 = 1;

% Define ODEs
odefun = @(t,y) [k1*X0 - k2*y(1); k3*X0/(1 + (y(1)/k4)^n) - k5*y(2)];

% Solve ODEs
[t_out, y] = ode45(odefun, t_sim, [0; 0]);

% Interpolate to match the size of avg_I1FFL
y_sim_I1FFL = interp1(t_out, y(:,2), t_sim);
end