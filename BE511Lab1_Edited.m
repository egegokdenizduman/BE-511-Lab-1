clearvars
clc
close all;
if ~exist('Figures2', 'dir')
    mkdir('Figures2');
end
set(groot, 'DefaultFigureWindowStyle', 'normal'); % new figures start undocked so the size sticks
voltage = linspace(2, 12, 6);
power = 0.05;
Vs = 10;
Is = linspace(0., 0.025, 6);
fontSize = 14; % pt at the 6.5 in width the figures go into the report
for i = 1:6
    current(i) = power/voltage(i);
    resistance(i) = voltage(i) / current(i);
    Rs(i) = Vs/current(i);
end

for i = 1:6
    voltage1(i) = -1*Is(i)*Rs(1) + Vs;
    voltage2(i) = -1*Is(i)*Rs(2) + Vs;
    voltage3(i) = -1*Is(i)*Rs(3) + Vs;
    voltage4(i) = -1*Is(i)*Rs(4) + Vs;
    voltage5(i) = -1*Is(i)*Rs(5) + Vs;
    voltage6(i) = -1*Is(i)*Rs(6) + Vs;
end

vSmooth = linspace(2, 12, 200); % fine sweep so the max power line is a smooth curve

figure(1);
clf;
set(gcf, 'WindowStyle', 'normal', 'Units', 'inches', 'Position', [1 1 6.5 5.75]); % taller to fit legend below
hold on;
plot(power ./ vSmooth, vSmooth, 'k-', 'LineWidth', 2.5);
plot(current, voltage, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7); % the 2 V steps
plot(Is, voltage1, '--', 'LineWidth', 2); % CHANGED: new line, 400 ohm is unsafe
plot(Is, voltage2, 'LineWidth', 2);
plot(Is, voltage3, 'LineWidth', 2);
plot(Is, voltage4, 'LineWidth', 2);
plot(Is, voltage5, 'LineWidth', 2);
plot(Is, voltage6, 'LineWidth', 2);
xlabel("Current (A)");
ylabel("Voltage (V)");
title("Maximum Power Line and Load Lines (V_s = 10 V)");
xlim([0, .025]);
ylim([0, 12]);
legend("Maximum Power Line (50 mW)", "Max Power at 2 V Steps", "R_s = 400 \Omega (unsafe)", ...
    "R_s = 800 \Omega", "R_s = 1200 \Omega", "R_s = 1600 \Omega", "R_s = 2000 \Omega", ...
    "R_s = 2400 \Omega", 'Location', 'southoutside', 'NumColumns', 2);
set(findall(gcf, '-property', 'FontWeight'), 'FontSize', fontSize, 'FontWeight', 'bold');
exportgraphics(gcf, 'Figures2/partA_load_lines.png', 'Resolution', 300);

%%

% part i
ThermometerTemp = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70];
ThermistorTemp = [25.62, 30.4, 36.37, 41.46, 44.48, 52.16, 57.4, 61.78, 66.92, 72.73];
ThermistorVoltage = [8.25, 8.05, 7.69, 7.36, 7.13, 6.61, 6.23, 5.92, 5.52, 5.14];
B = 3380;
T0 = 25 + 273.15;
R0 = 10000;
RsMeasured = 1968.6; % measured source resistance used in Part D

% expected curves from the beta equation and voltage divider
tSmooth = linspace(25, 70, 200);
RSmooth = R0 * exp(B ./ (tSmooth + 273.15) - B / T0);
VSmooth = Vs * RSmooth ./ (RSmooth + RsMeasured);

figure(2);
clf;
set(gcf, 'WindowStyle', 'normal', 'Units', 'inches', 'Position', [1 1 6.5 4.5]); % undock so the size sticks
hold on;
plot(ThermometerTemp, ThermistorVoltage, 'o', 'MarkerFaceColor', [0 0.447 0.741], 'MarkerSize', 8);
plot(tSmooth, VSmooth, '-', 'LineWidth', 2);
xlim([20, 75]);
xlabel('Thermometer Temperature (\circC)');
ylabel('Thermistor Voltage (V)');
title('Thermistor Voltage vs. Thermometer Temperature');
legend('Observed Voltage', 'Expected Voltage', 'Location', 'northeast');
set(findall(gcf, '-property', 'FontWeight'), 'FontSize', fontSize, 'FontWeight', 'bold');
exportgraphics(gcf, 'Figures2/partD_voltage_vs_temp.png', 'Resolution', 300);

% part ii
for i = 1:length(ThermistorTemp)
    exponentTerm = (B/(ThermistorTemp(i)+273.15)) - B/T0;
    ThermistorResistance(i) = R0 * exp(exponentTerm);
    exponentTerm2 = (B/(ThermometerTemp(i)+273.15)) - B/T0;
    ExpectedThermistorResistance(i) = R0 * exp(exponentTerm2);
end

figure(3);
clf;
set(gcf, 'WindowStyle', 'normal', 'Units', 'inches', 'Position', [1 1 6.5 4.5]); % undock so the size sticks
hold on;
plot(ThermometerTemp, ThermistorResistance, 'o', 'MarkerFaceColor', [0 0.447 0.741], 'MarkerSize', 8);
plot(tSmooth, RSmooth, '-', 'LineWidth', 2); % smooth expected curve
xlim([20, 75]);
xlabel('Thermometer Temperature (\circC)');
ylabel('Thermistor Resistance (\Omega)');
title('Thermistor Resistance vs. Thermometer Temperature');
legend('Observed Resistance', 'Expected Resistance', 'Location', 'northeast');
set(findall(gcf, '-property', 'FontWeight'), 'FontSize', fontSize, 'FontWeight', 'bold');
exportgraphics(gcf, 'Figures2/partD_resistance_vs_temp.png', 'Resolution', 300);

% part iii
figure(4);
clf;
set(gcf, 'WindowStyle', 'normal', 'Units', 'inches', 'Position', [1 1 6.5 4.5]); % undock so the size sticks
hold on;
plot(ThermometerTemp, ThermistorTemp, 'o', 'MarkerFaceColor', [0 0.447 0.741], 'MarkerSize', 8);
plot(ThermometerTemp, ThermometerTemp, 'k--', 'LineWidth', 2);
xlim([20, 75]);
xlabel('Thermometer Temperature (\circC)');
ylabel('Thermistor Temperature (\circC)');
title('Thermistor Temperature vs. Thermometer Temperature');
legend('Thermistor Temperature', 'Ideal Thermistor Temperature', 'Location', 'northwest');
set(findall(gcf, '-property', 'FontWeight'), 'FontSize', fontSize, 'FontWeight', 'bold');
exportgraphics(gcf, 'Figures2/partD_thermistor_vs_thermometer.png', 'Resolution', 300);

% part iv
for i = 1:length(ThermistorTemp)
    residual(i) = ThermistorTemp(i) - ThermometerTemp(i);
end

figure(5);
clf;
set(gcf, 'WindowStyle', 'normal', 'Units', 'inches', 'Position', [1 1 6.5 4.5]); % undock so the size sticks
hold on;
stem(ThermometerTemp, residual, 'filled', 'LineWidth', 2, 'MarkerSize', 8); % line from 0 shows size of each error
xlim([20, 75]);
xlabel('Thermometer Temperature (\circC)');
ylabel('Residual (\circC)');
title('Residual (Thermistor - Thermometer)');
set(findall(gcf, '-property', 'FontWeight'), 'FontSize', fontSize, 'FontWeight', 'bold');
exportgraphics(gcf, 'Figures2/partD_residuals.png', 'Resolution', 300);
