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
salineControl = xlsread(Summarized_files(1).name);
lowNic = xlsread(Summarized_files(2).name);
highNic = xlsread(Summarized_files(3).name);

global random_order;
n = 58; % Set the maximum number
random_order = randperm(n)'; % Generate a random permutation

% Normalize and sort data
[salineControl, salineControl_sorted_oriData, salineControl_random_order] = normalizeAndSort(salineControl, 'descend');
[lowNic, lowNic_sorted_oriData, lowNic_random_order] = normalizeAndSort(lowNic, 'descend');
[highNic, highNic_sorted_oriData, highNic_random_order ] = normalizeAndSort(highNic, 'ascend');

% Initialize a new figure
figure;
% Set figure size
set(gcf, 'Position', [10, 10, 1200, 900]);

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
print(gcf, 'Heatmap.pdf', '-dpdf', '-fillpage');

figure
% Overlay all three groups in one plot
% subplot(3,3,7)
% plotOverlay(salineControl, lowNic, highNic);

hold on;

% Function to plot mean and SEM
plotMeanAndSEM(salineControl, [0, 0, 0]); % Black for Saline control
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
[duration_increase_lowNic, duration_decrease_lowNic] = calcDuration(lowNic, 300);
[duration_increase_highNic, duration_decrease_highNic] = calcDuration(highNic, 300);

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


% Scatter Plot with Line of Best Fit
figure;

scatter(common_durations_lowNic, common_durations_highNic, 'k');
axis square;
xlim([200 900]); % Set x-axis limits
ylim([200 900]); % Set y-axis limits

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

%%
% New figure for 2D plot of High Nicotine group with both leftward and upward shifts for each trial
figure;
% Set figure size
set(gcf, 'Position', [10, 10, 1500, 600]);
subplot(1,3,1)
hold on;
% Plot 2D raw traces for High Nic with smoothing, time range from 150s to 450s, and both types of shifts
plotRawTraces2DWithBothShiftsAndLabels(salineControl, [0, 0, 0], 150, 450);
xlabel('Time (s)');
ylabel('Normalized Firing Rate');
title('Saline Control');

subplot(1,3,2)
hold on;
% Plot 2D raw traces for High Nic with smoothing, time range from 150s to 450s, and both types of shifts
plotRawTraces2DWithBothShiftsAndLabels(lowNic, [1, 0.5, 0], 150, 450);
xlabel('Time (s)');
ylabel('Normalized Firing Rate');
title('Low Nic');

subplot(1,3,3)
hold on;
% Plot 2D raw traces for High Nic with smoothing, time range from 150s to 450s, and both types of shifts
plotRawTraces2DWithBothShiftsAndLabels(highNic, [0.5, 0, 0.5], 150, 450);
xlabel('Time (s)');
ylabel('Normalized Firing Rate');
title('High Nic');

% grid on;


% Save figure as PDF
saveas(gcf, 'HighNic2DShiftedBothWaysWithLabels.pdf', 'pdf');

% Function to plot 2D raw traces with time range, smoothing, both types of shifts, and trial labels
function plotRawTraces2DWithBothShiftsAndLabels(data, color, startTime, endTime)
% Assuming 1 sample per second, adjust this as per your actual sampling rate
samplingRate = 1; % Change this to match your actual data's sampling rate
startIndex = startTime * samplingRate + 1;
endIndex = endTime * samplingRate + 1;

% Smoothing the data using a moving average
windowSize = 5; % Change this value to adjust the smoothness
smoothedData = movmean(data, windowSize, 2);

timeVector = linspace(startTime, endTime, endIndex - startIndex + 1);
verticalShiftAmount = 0.15; % Change this value to adjust vertical spacing between trials
horizontalShiftAmount = 2; % Change this value to adjust horizontal spacing between trials

for i = 1:size(data, 1)
    % Applying both vertical and horizontal shifts
    shiftedTime = timeVector + (i-1) * horizontalShiftAmount;
    shiftedData = smoothedData(i, startIndex:endIndex) + (i-1) * verticalShiftAmount;
    plot(shiftedTime, shiftedData, 'Color', color);

    % Adding text labels for each trial number
    trialNumberText = sprintf('Trial %d', i);
    % text(shiftedTime(end) + 1, shiftedData(end), trialNumberText, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 8, 'Color', color);
end
axis square;
ylim([0 10]);
end

hold off;





% Helper function to normalize and sort data
function [sortedData, sorted_oriData, sorted_Data_random_order] = normalizeAndSort(data, sortDirection)
global random_order
% Normalize each row to its mean from 1 to 200
Ori_data = data;
for i = 1:size(data, 1)
    meanValue = mean(data(i, 1:200));
    Mean_data(i, :) = Ori_data(i, :) / meanValue;
end

% Sort rows based on average from bins 350 to end
[~, sortIdx] = sort(mean(Mean_data(:, 350:end), 2), sortDirection);
sortedData = Mean_data(sortIdx, :);
sorted_oriData = Ori_data (sortIdx, :);
sorted_Data_random_order = sorted_oriData(random_order,:);
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
meanFirst300 = mean(data(:, 1:300), 2);
mean300to600 = mean(data(:, 300:1200), 2);
unitsWithReduction = mean300to600 < 0.8 * meanFirst300;
unitsWithIncreasing = mean300to600 > 1.2 * meanFirst300;
reductionRatio = sum(unitsWithReduction) / size(data, 1);
increaseRatio = sum(unitsWithIncreasing) / size(data, 1);

increased = increaseRatio * 100;
decreased = reductionRatio * 100;
unchanged = 100 - increased - decreased;

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
% Adjust the timeVector to match the new data dimensions
timeVector = linspace(0, length(meanData)-1, length(meanData));

% Ensure the vectors for the fill function have matching dimensions
upperSEM = meanData + SEM;
lowerSEM = meanData - SEM;

fill([timeVector, fliplr(timeVector)], [lowerSEM, fliplr(upperSEM)], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(timeVector, meanData, 'Color', color, 'LineWidth', 2); % Plot the mean on top

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
