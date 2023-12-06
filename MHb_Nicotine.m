% Clear workspace, close all figures, and clear command window
clc;
close all;
clear all;

% List all xlsx files in the current directory
Summarized_files = dir('*.xlsx');

% Initialize variables to store data
salineControl = [];
lowNic = [];
highNic = [];

% Load data from 'Original_freq' spreadsheet for Saline control, Low Nic, and High Nic
salineControl = xlsread(Summarized_files(1).name, 'Original_freq');
lowNic = xlsread(Summarized_files(2).name, 'Original_freq');
highNic = xlsread(Summarized_files(3).name, 'Original_freq');

% Normalize and sort data
salineControl = normalizeAndSort(salineControl, 'descend');
lowNic = normalizeAndSort(lowNic, 'descend');
highNic = normalizeAndSort(highNic, 'ascend');

% ... [You can continue with your plotting code here]

% ... [Previous code to load and normalize data]

% Initialize a new figure
figure;
% Set figure size
set(gcf, 'Position', [100, 100, 1600, 1200]);

% Plotting heatmaps
subplot(2,3,1)
plotHeatmap(salineControl, 'Saline Control');
subplot(2,3,2)
plotHeatmap(lowNic, 'Low Nic');
subplot(2,3,3)
plotHeatmap(highNic, 'High Nic');

% Plotting pie charts
subplot(2,3,4)
plotPieChart(salineControl, 'Saline Control');
subplot(2,3,5)
plotPieChart(lowNic, 'Low Nic');
subplot(2,3,6)
plotPieChart(highNic, 'High Nic');

% Save figure as PDF (or you can choose EPS)
saveas(gcf, 'Heatmap.pdf', 'pdf');

figure
% Overlay all three groups in one plot
% subplot(3,3,7)
% plotOverlay(salineControl, lowNic, highNic);

hold on;

% Function to plot mean and SEM
plotMeanAndSEM(salineControl, 'k'); % Black for Saline control
plotMeanAndSEM(lowNic, [1, 0.5, 0]); % Orange for Low Nic
plotMeanAndSEM(highNic, [0.5, 0, 0.5]); % Purple for High Nic

ylim([0.5,1.5])
axis square;
xlabel('Time (s)');
ylabel('Normalized Firing Rate');
title('Overlay of All Groups');
legend('Saline control', 'Low Nic', 'High Nic', 'Location', 'best');
hold off;

% Save figure as PDF (or you can choose EPS)
saveas(gcf, 'MyFigure.pdf', 'pdf');



% Calculate duration of firing change for each condition and each unit
[duration_increase_lowNic, duration_decrease_lowNic] = calcDuration(lowNic, 30);
[duration_increase_highNic, duration_decrease_highNic] = calcDuration(highNic, 30);

Thresthold_1 = 20;

% Find units that show an increase in low Nic and decrease in high Nic
increasing_units_lowNic = find(duration_increase_lowNic > Thresthold_1);
decreasing_units_highNic = find(duration_decrease_highNic > Thresthold_1);

common_units = intersect(increasing_units_lowNic, decreasing_units_highNic);

% Pearson correlation between the duration of increase in low Nic and duration of decrease in high Nic for common units
common_durations_lowNic = duration_increase_lowNic(common_units);
common_durations_highNic = duration_decrease_highNic(common_units);

[R, P] = corrcoef(common_durations_lowNic, common_durations_highNic); % Pearson
% [R, P] = corr(common_durations_lowNic, common_durations_highNic, 'Type', 'Spearman'); % Spearman


% % Scatter Plot with Mean and SEM
% figure;
% scatter(common_durations_lowNic, common_durations_highNic, 'filled');
% 
% axis square;
% hold on;
% 
% % Calculate the mean and SEM
% mean_lowNic = mean(common_durations_lowNic);
% mean_highNic = mean(common_durations_highNic);
% 
% sem_lowNic = std(common_durations_lowNic) / sqrt(length(common_durations_lowNic));
% sem_highNic = std(common_durations_highNic) / sqrt(length(common_durations_highNic));
% 
% % Plot mean + SEM lines
% line([mean_lowNic, mean_lowNic], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% line(xlim, [mean_highNic, mean_highNic], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% 
% % Plot mean + SEM points
% errorbar(mean_lowNic, mean_highNic, sem_highNic, 'horizontal', 'r', 'LineWidth', 2);
% errorbar(mean_lowNic, mean_highNic, sem_lowNic, 'vertical', 'r', 'LineWidth', 2);
% 
% xlabel('Duration of Increase in Low Nicotine (s)');
% ylabel('Duration of Decrease in High Nicotine (s)');
% title(['Pearson Correlation: R = ' num2str(R(1,2)) ', P = ' num2str(P(1,2))]);
% grid on;
% hold off;
% 


% Scatter Plot with Line of Best Fit
figure;

scatter(common_durations_lowNic, common_durations_highNic, 'k');
axis square;
xlim([20 90]); % Set x-axis limits
ylim([20 90]); % Set y-axis limits

xlabel('Duration of Increase in Low Nicotine (s)');
ylabel('Duration of Decrease in High Nicotine (s)');
title(['Pearson Correlation: R = ' num2str(R(1,2)) ', P = ' num2str(P(1,2))]);

% Calculate and plot the line of best fit
coeffs = polyfit(common_durations_lowNic, common_durations_highNic, 1);
fitLine = polyval(coeffs, common_durations_lowNic);
hold on; % keep the current scatter plot
plot(common_durations_lowNic, fitLine, 'r-', 'LineWidth', 2); % plot the line of best fit
hold off;
% Save figure as PDF (or you can choose EPS)
saveas(gcf, 'Pearson Correlation.pdf', 'pdf');

% legend('Data Points', 'Line of Best Fit', 'Location', 'best');

% grid on;







% Helper function to normalize and sort data
function [sortedData] = normalizeAndSort(data, sortDirection)
    % Normalize each row to its mean from 1 to 20
    for i = 1:size(data, 1)
        meanValue = mean(data(i, 1:20));
        data(i, :) = data(i, :) / meanValue;
    end
    
    % Sort rows based on average from bins 35 to end
    [~, sortIdx] = sort(mean(data(:, 35:end), 2), sortDirection);
    sortedData = data(sortIdx, :);
end

% Helper function to plot heatmap
function plotHeatmap(data, titleStr)
    colormap_clim = [0 5];
    imagesc(data, colormap_clim); 
    axis square;
    colormap(gca,'hot');
    colorbar;
    xlabel('Time');
    ylabel('Units');
    title(['Heatmap of ' titleStr]);
end

% Helper function to plot pie chart
function plotPieChart(data, titleStr)
    meanFirst30 = mean(data(:, 1:30), 2);
    mean30to60 = mean(data(:, 30:120), 2);
    unitsWithReduction = mean30to60 < 0.8 * meanFirst30;
    unitsWithIncreasing = mean30to60 > 1.2 * meanFirst30;
    reductionRatio = sum(unitsWithReduction) / size(data, 1);
    increaseRatio = sum(unitsWithIncreasing) / size(data, 1);
    
    increased = increaseRatio * 100
    decreased = reductionRatio * 100
    unchanged = 100 - increased - decreased
    
    labels = {'Increased', 'Decreased', 'Unchanged'};
    custom_colormap = [1, 0.5, 0.5; 0.2, 0.2, 0.2; 1, 1, 1];
    pie([increased, decreased, unchanged], labels);
    colormap(gca, custom_colormap);
    title(['Pie Chart of ' titleStr]);
end

% Helper function to overlay plots
function plotOverlay(salineControl, lowNic, highNic)
    hold on;
    plot(mean(salineControl, 1), 'Color', 'k');
    plot(mean(lowNic, 1), 'Color', [0.6, 1, 0.6]);
    plot(mean(highNic, 1), 'Color', [1, 0.6, 0.6]);
    hold off;
    
    legend('Saline Control', 'Low Nic', 'High Nic');
    xlabel('Time');
    ylabel('Mean Firing Rate');
    title('Overlay of All Groups');
end

% Function to plot mean and SEM
function plotMeanAndSEM(data, color)
    meanData = mean(data, 1);
    SEM = std(data, 0, 1) / sqrt(size(data, 1));
    timeVector = 0:10:1190; % Assuming your data has 120 bins
    fill([timeVector, fliplr(timeVector)], [meanData - SEM, fliplr(meanData + SEM)], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(timeVector, meanData, 'Color', color, 'LineWidth', 2);
end




% Helper function to calculate the duration of firing increase and decrease
function [duration_increase, duration_decrease] = calcDuration(data, baselineTime)
    nUnits = size(data, 1);
    duration_increase = zeros(nUnits, 1);
    duration_decrease = zeros(nUnits, 1);

    for i = 1:nUnits
        baseline = mean(data(i, 1:baselineTime));
        post_data = data(i, baselineTime+1:end);
        Threshold = 0.2;
        duration_increase(i) = sum(post_data > baseline + Threshold);
        duration_decrease(i) = sum(post_data < baseline - Threshold );
    end
end
