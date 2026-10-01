clearvars
clc
%% Normalization and data procesisng
% Constants
Vs = 10;
B = 3380;
T0 = 25 + 273.15;
R0 = 10000;
Rs = [1968.6, 3966.4, 994.3]; % measured source resistance for each setup
threshold = 0.1;              % deg C per sample to locate sumbersion
dt = 0.110505;                % sample period in s
fontSize = 14;                % pt at the 6.5 in width the figures go into the report
dataFolder = 'Data';          % folder with the exported test files

% load data and find the dunk and removal sample of each test
for p = 1:3
    for t = 1:3
        data = readmatrix(fullfile(dataFolder, sprintf('voltage part %d test %d.txt', p, t)), 'NumHeaderLines', 1);
        voltage{p, t} = data(:, 2);
        resistance = voltage{p, t} * Rs(p) ./ (Vs - voltage{p, t});
        temp{p, t} = 1 ./ (log(resistance / R0) / B + 1 / T0) - 273.15;

        % dunk
        dT = diff(temp{p, t});
        [~, s] = max(dT);
        while dT(s - 1) > threshold
            s = s - 1;
        end

        % removal
        [~, e] = min(dT);
        while dT(e - 1) < -threshold
            e = e - 1;
        end

        dunk(p, t) = s;
        removal(p, t) = e;
    end
end

%% Plot

for p = 1:3
    % samples kept before, in and out of water
    before = min(dunk(p, :)) - 1;
    plateau = min(removal(p, :) - dunk(p, :));
    after = inf;
    for t = 1:3
        after = min(after, length(temp{p, t}) - removal(p, t) + 1);
    end
    x = (-before : plateau + after - 1)' * dt; % time from dunk
    cutTime = plateau * dt; % where the plateau was cut

    tempAligned = zeros(length(x), 3);
    voltAligned = zeros(length(x), 3);
    for t = 1:3
        s = dunk(p, t);
        e = removal(p, t);
        idx = [s - before : s + plateau - 1, e : e + after - 1];
        tempAligned(:, t) = temp{p, t}(idx);
        voltAligned(:, t) = voltage{p, t}(idx);
    end

% Plot fr   
    tempColor = [0 0.447 0.741];    % matlab blue
    voltColor = [0.85 0.325 0.098]; % matlab red
    figure(p);
    clf; % debug thingy
    set(gcf, 'Units', 'inches', 'Position', [1 1 6.5 7.5]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');

    nexttile;
    yyaxis left;
    plotMeanStd(x, tempAligned, tempColor);
    ylabel('Temperature (\circC)');
    yyaxis right;
    plotMeanStd(x, voltAligned, voltColor);
    ylabel('Voltage (V)');
    xline(0, 'k--', 'LineWidth', 1.5);
    xline(cutTime, 'k:', 'LineWidth', 1.5);
    xlabel('Time from Immersion (s)');
    title(sprintf('Temperature and Voltage (R_s = %.1f \\Omega)', Rs(p)));
    xlim([x(1), x(end)]);

    nexttile;
    yyaxis left;
    hTemp = plotMeanStd(x, tempAligned - tempAligned(1, :), tempColor);
    ylabel('\DeltaTemperature (\circC)');
    yyaxis right;
    hVolt = plotMeanStd(x, voltAligned - voltAligned(1, :), voltColor);
    ylabel('\DeltaVoltage (V)');
    hDunk = xline(0, 'k--', 'LineWidth', 1.5);
    hCut = xline(cutTime, 'k:', 'LineWidth', 1.5);
    xlabel('Time from Immersion (s)');
    title(sprintf('Change from First Sample (R_s = %.1f \\Omega)', Rs(p)));
    xlim([x(1), x(end)]);
    legend([hTemp, hVolt, hDunk, hCut], 'Temperature \pm 1 SD', 'Temperature Mean', ...
        'Voltage \pm 1 SD', 'Voltage Mean', 'Immersion', 'Removal', ...
        'Location', 'southoutside', 'NumColumns', 3); % one legend for both panels, outside the data

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);
    exportgraphics(gcf, sprintf('Figures2/partE_mean_setup_%d.png', p), 'Resolution', 300);

end

%% Tau

for p = 1:3
    figure(p + 3);
    clf; % debbug again
    set(gcf, 'Units', 'inches', 'Position', [1 1 6.5 7.5]);
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact');
    for t = 1:3
        T = temp{p, t};
        tm = ((1:length(T))' - dunk(p, t)) * dt; % time from dunk
        s = dunk(p, t);
        e = removal(p, t);

        initialValue = T(s);
        finalValue = mean(T(e - 9:e)); % water plateau before removal
        target = initialValue + (finalValue - initialValue) * (1 - 1/exp(1));

        % interpolation
        c = s;
        while T(c) < target
            c = c + 1;
        end
        endTime = interp1([T(c - 1), T(c)], [tm(c - 1), tm(c)], target);
        tau(p, t) = endTime - tm(s);

        % Plot
        nexttile;
        hold on;
        plot(tm, T, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        plot(tm(s), T(s), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
        plot(endTime, target, 'rx', 'MarkerSize', 14, 'LineWidth', 2.5);
        yline(target, 'r:', 'LineWidth', 1.5);
        xlim([-1, 5]);
        xlabel('Time from Immersion (s)');
        ylabel('Temperature (\circC)');
        title(sprintf('Test %d, \\tau = %.3f s', t, tau(p, t)));
        if t == 1
            legend('Temperature', 'Start', 'End (63.2%)', '63.2% Level', 'Location', 'southeast');
        end
    end
    title(tl, sprintf('R_s = %.1f \\Omega', Rs(p)));
    fprintf('Rs = %.1f ohms: tau = %.3f +/- %.3f s\n', Rs(p), mean(tau(p, :)), std(tau(p, :)));
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);
    exportgraphics(gcf, sprintf('Figures2/partE_tau_dunk_setup_%d.png', p), 'Resolution', 300);
end

%% Voltage vs Temperature

figure(7);
clf; % debug thingy
set(gcf, 'Units', 'inches', 'Position', [1 1 6.5 4.5]);
hold on;
for p = 1:3
    vAll = vertcat(voltage{p, :});
    tAll = vertcat(temp{p, :});
    [vAll, i] = sort(vAll); % sort by voltage so the line is a curve
    plot(vAll, tAll(i), '-', 'LineWidth', 2); % all 3 tests
end
grid on;
xlabel('Thermistor Voltage (V)');
ylabel('Temperature (\circC)');
title('Temperature vs. Thermistor Voltage for Each R_s');
legend(sprintf('R_s = %.1f \\Omega', Rs(1)), sprintf('R_s = %.1f \\Omega', Rs(2)), ...
    sprintf('R_s = %.1f \\Omega', Rs(3)), 'Location', 'southwest');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);
exportgraphics(gcf, 'Figures2/partE_temp_vs_voltage.png', 'Resolution', 300);

% mean +/- std
function h = plotMeanStd(x, Y, color)
    mu = mean(Y, 2);
    sd = std(Y, 0, 2);
    hold on;
    hFill = fill([x; flipud(x)], [mu - sd; flipud(mu + sd)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hMean = plot(x, mu, '-', 'Color', color, 'LineWidth', 2);
    h = [hFill, hMean];
end
