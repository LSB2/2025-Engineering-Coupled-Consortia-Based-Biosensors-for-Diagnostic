% Single Input
clc; close all; clear all;

UI1 = 0.003;
UI2 = 0.004;
N01 = 1;
N02 = 1;

NI0 = 0.1;
tspan = 0:1:2400;

init_NI1 = NI0; % initial condition vector
odes = @(t_NI1, NI1) UI1*NI1*(1-NI1/N02);
[t_NI1, NI1] = ode45(odes, tspan, init_NI1);

init_NI2 = NI0;
start_time = 0 * 60;  % start time is set at 0 for simplicity

% No fluctuations for NI2
odes_without_fluctuation = @(t_NI2, NI2) (t_NI2 >= start_time) * (UI2 * NI2 * (1 - NI2 / N02));
[t_NI2_without_fluctuation, NI2_without_fluctuation] = ode45(odes_without_fluctuation, tspan, init_NI2);

% Introduce noise to NI2
noise_amplitude = 0.05;
%noise_amplitude = 0.003;
odes_with_noise = @(t_NI2, NI2) (t_NI2 >= start_time) * (UI2 * NI2 * (1 - NI2 / N02) + noise_amplitude * randn);
[t_NI2_with_noise, NI2_with_noise] = ode45(odes_with_noise, tspan, init_NI2);

%Figure1 - Population for both strains
fig = figure;
%fig.Position = [55 550 560 420];
fig.Position = [55 550 500 375];
hold on
plot(t_NI1./60, NI1, 'LineWidth', 2)
plot(t_NI2_without_fluctuation./60, NI2_without_fluctuation, '--', 'LineWidth', 2) % dashed line for no fluctuations
plot(t_NI2_with_noise./60, NI2_with_noise, ':', 'LineWidth', 2) % dotted line for with noise
legend('Input 1','Input 2 without fluctuations','Input 2 with noise', 'Location', 'southeast')
xlabel('Cell Generations')
ylabel('OD_{600}')
title('Population of the Input Strain')
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
box on;
hold off

% Parameters
a_In1 = 1;
a_In2 = 1;
K_In1 = 1;
K_In2 = 1;
b_In1 = 0.05;
b_In2 = 0.15;

X0 = 1; % concentration of X
Y0 = 0; % initial concentration of Y
Z0 = 0; % initial concentration of Z

tspan = 0:40/2400:40; % time span for the simulation
init = [Y0; Z0]; % initial condition vector

% Initial conditions
In10 = 0.1; % initial concentration of
In20 = 0.1;

odes_In1 = @(t, In1) a_In1*X0./(K_In1+X0) - b_In1*In1;
[t_In1, In1] = ode45(odes_In1, tspan, Z0);
In1 = (In1 - min(In1)) / (max(In1) - min(In1)); % Normalize

odes_In2 = @(t, In2) a_In2*X0./(K_In2+X0) - b_In2*In2;
[t_In2, In2] = ode45(odes_In2, tspan, Z0);
In2 = (In2 - min(In2)) / (max(In2) - min(In2)); % Normalize

%Figure – Input Signal f(I1) and f(I2)
fig = figure;
%fig.Position = [630 550 560 420];
fig.Position = [55 550 500 375];
hold on
plot(t_In1, In1, 'LineWidth', 2)
plot(t_In2, In2, 'LineWidth', 2)
legend('f(I_1)','f(I_2)', 'Location', 'best')
%title('Input Signal')
xlabel('Cell Generations')
ylabel('Input Signal')
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
box on
hold off

fin1_total = (NI1.*In1 + NI2_without_fluctuation.*In2);
max_fin1 = max(fin1_total);
Norm_fin1_total = fin1_total / max_fin1;

fin2_total = (NI1.*In1+NI2_with_noise.*In2);
max_fin2 = max(fin2_total);
Norm_fin2_total = fin2_total / max_fin2;


fig = figure;
fig.Position = [55 40 2*500 375];
subplot(1,2,1)
hold on
plot(t_NI1./60, Norm_fin1_total, 'LineWidth', 2)
plot(t_NI1./60, Norm_fin2_total, 'LineWidth', 2)
legend('Without fluctuation','With fluctuation', 'Location', 'northwest')
xlabel('Cell Generations')
title('Independent Case')
ylabel('Normalized Output')
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
box on;
hold off


f1_total = (NI1.*In1+NI2_without_fluctuation.*In2)./(NI1+NI2_without_fluctuation);
max_f1 = max(f1_total);
Norm_f1_total = f1_total / max_f1;

f2_total = (NI1.*In1+NI2_with_noise.*In2)./(NI1+NI2_with_noise);
max_f2 = max(f2_total);
Norm_f2_total = f2_total / max_f2;

subplot(1,2,2)
hold on
plot(t_NI1./60, Norm_f1_total, 'LineWidth', 2)
plot(t_NI1./60, Norm_f2_total, 'LineWidth', 2)
legend('Without fluctuation','With fluctuation', 'Location', 'southeast')
xlabel('Cell Generations')
title('Normalized Output')
title('Shared Signal Case')
set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
box on;
hold off


St=[10 1 0.1];
n=length(St);

fig = figure;
fig.Position = [55 40 3*500 375];
colors = lines(n); % Using MATLAB's lines colormap to generate distinct colors
legend_entries = cell(1, 2*n); % Initializing for storing legend info

for i=1:n
    subplot(1,n,i)
    fsh1_total = St(i).*(NI1.*In1 + NI2_without_fluctuation.*In2) ./ (St(i) + NI1 + NI2_without_fluctuation);
    Norm_fsh1_total = fsh1_total / max(fsh1_total);

    fsh2_total = St(i).* (NI1.*In1 + NI2_with_noise.*In2) ./ (St(i) + NI1 + NI2_with_noise);
    Norm_fsh2_total = fsh2_total / max(fsh2_total);
    hold on
    plot(t_NI1./60, Norm_fsh1_total, 'LineWidth', 2, 'Color', colors(i, :))
    plot(t_NI1./60, Norm_fsh2_total, '--', 'LineWidth', 2, 'Color', colors(i, :)) % Using dashed line for the second plot

    legend_entries = cell(1, 2); % Local legend entries for this subplot
    legend_entries{1} = sprintf('without fluctuation at St = %g', St(i));
    legend_entries{2} = sprintf('with fluctuation at St = %g', St(i));
    legend(legend_entries, 'Location', 'southeast'); % Show legend to differentiate curves

    hold off;
    xlabel('Cell Generations')
    ylabel('Normalized Output')
    %title('Shared Signal Case')
    set(gca, 'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    box on
end

