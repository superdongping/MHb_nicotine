load('Single_unit_58neurons_ori_firing_data.mat')
% Define the baseline period (first 50 seconds)
baselinePeriod = 50;  % in seconds
baselineBins = 1:baselinePeriod;  % Corresponding bins for the baseline period

% Initialize an array to store the reduction in firing rates for each neuron
numNeurons = 58;
firingRateReductions = zeros(numNeurons, 1);

% Calculate the reduction in firing rate for each neuron
for i = 1:numNeurons
    % Calculate the baseline firing rate (mean firing rate during the baseline period)
    baselineFiringRate = mean(Neuron_FiringRates(i, baselineBins));
    
    % Calculate the post-baseline mean firing rate
    postBaselineFiringRate = mean(Neuron_FiringRates(i, baselinePeriod+1:end));
    
    % Calculate the reduction as a percentage of the baseline firing rate
    firingRateReductions(i) = (baselineFiringRate - postBaselineFiringRate) / baselineFiringRate;
end

% Set a 20% reduction threshold and find neurons that meet the criteria
reductionThreshold = 0.2;
selectedNeurons = find(firingRateReductions > reductionThreshold);

% Sort the selected neurons based on the magnitude of reduction
[~, sortedIndices] = sort(firingRateReductions(selectedNeurons), 'descend');
sortedNeurons = selectedNeurons(sortedIndices);

% Select a few neurons from the sorted list (choose the best ones based on your judgment)
selectedBestNeurons = sortedNeurons([1, 2, 5,6 ,7,9,11,12,13,16,18,21,23,25]);  % You can manually change these indices

% Create a new figure to plot the selected best neurons
figure('Units', 'normalized', 'Position', [0.1 0.2 0.8 0.5]);  % Adjust the figure size for better view
hold on;

% Define a colormap for firing rates (e.g., jet colormap)
cmap = jet(100);  % Creates a colormap with 100 color levels

% Get the minimum and maximum firing rates across all neurons for scaling the colormap
minFiringRate = min(Neuron_FiringRates(:));  % Minimum firing rate
maxFiringRate = max(Neuron_FiringRates(:));  % Maximum firing rate

% Ensure that the color axis has proper limits
caxis([minFiringRate maxFiringRate]);

% Plot the selected neurons
for n = 1:length(selectedBestNeurons)
    i = selectedBestNeurons(n);  % Get the neuron index from sorted list
    spks = Neuron_TimeStamp_sum{i};  % Get spike times for this neuron
    
    for j = 1:length(spks)
        % Determine the corresponding bin for the spike
        spikeTime = spks(j);
        bin = floor(spikeTime) + 1;  % Convert spike time to bin index

        % Get the firing rate for this neuron and bin
        firingRate = Neuron_FiringRates(i, bin);

        % Normalize the firing rate to a value between 1 and 100 based on real firing rate
        normRate = min(max(round((firingRate - minFiringRate) / (maxFiringRate - minFiringRate) * 99) + 1, 1), 100); 

        % Select the corresponding color from the colormap
        color = cmap(normRate, :);

        % Plot spike with the appropriate color based on real firing rate
        plot([spikeTime spikeTime], [n-1+0.1 n-0.1], 'Color', color, 'LineWidth', 0.1);
    end
end

xlabel('Time (s)');
ylabel('Neuron Number');
title('Raster Plot');
xlim([ts tf]);
ylim([0 length(selectedBestNeurons)]);

% Add a colorbar to indicate the firing rate color scale
colormap(cmap);

% Make sure the colorbar is displayed
colorbar('eastoutside');  % Force the colorbar to the outside of the plot

% Set the color axis based on the actual firing rate range
caxis([0 maxFiringRate-5]);

% Label the colorbar with real firing rates
cb = colorbar;
cb.Label.String = 'Firing Rate (Hz)';

% Save the plot as a PDF file that preserves vector graphics
print('SelectedNeuronsPlot', '-dpdf', '-painters', '-r1200');



