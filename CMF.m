clc;
clear;
close all;

% ------------------------- Loead DATA
% Load the HDF5 file
  filePath = 'Taftan_timeseries_Corrected.h5';


% Load the 'timeseries' dataset
timeseries_data = h5read(filePath, '/timeseries');

% Load metadata for coordinate calculations from the HDF5 file
x_first = str2double(h5readatt(filePath, '/', 'X_FIRST'));
x_step = str2double(h5readatt(filePath, '/', 'X_STEP'));
y_first = str2double(h5readatt(filePath, '/', 'Y_FIRST'));
y_step = str2double(h5readatt(filePath, '/', 'Y_STEP'));
ncols = size(timeseries_data, 2); % Assuming data is [nrows x ncols x ntimes]
nrows = size(timeseries_data, 1);

% -------------------- Select POINT (e.g Ref. point)
% % Get the time series for the specific point  (REF Point)
% target_lat = 28.62662;
% target_lon = 61.14602;
% 
% target_lat =  28.601580;   % summit   p1
% target_lon =  61.134127;

% 
% target_lat = 28.600179;   % west    p2
% target_lon =  61.128392;
% 
% target_lat =  28.601775;    % East     p3
% target_lon =  61.138594;

% 
% target_lat =  28.599648;    % far
% target_lon =   61.121271;

%target_lat =   28.599183;    % west_second
%target_lon = 61.130367;

% % % % P5
% target_lat = 28.606028     % p4
% target_lon = 61.126624;

% % P6
target_lat = 28.590772;        % p5
target_lon = 61.128804;

% Calculate indices for the target point
target_x_idx = round((target_lon - x_first) / x_step);
target_y_idx = round((target_lat - y_first) / y_step);

% % Ensure indices are within valid range
% target_x_idx = max(1, min(size(timeseries_data, 1), target_x_idx));
% target_y_idx = max(1, min(size(timeseries_data, 2), target_y_idx));

% Extract the time series for the specific point
point_time_series = timeseries_data(target_x_idx, target_y_idx, :);         % (Row(X) , Col(Y)) ?! (664, 469)
point_time_series_cm = squeeze(point_time_series) * 100; % Convert to cm


% Define the bounding box (lower left and upper right corners)
% lower_left_lon = 61.08;   
% upper_right_lon = 61.18;
% lower_left_lat = 28.55;
% upper_right_lat = 28.66;

lower_left_lon = 61.11;  
lower_left_lat = 28.58;
upper_right_lon = 61.15;
upper_right_lat = 28.62;

% Calculate indices for the bounding box
x_idx_min = round((lower_left_lon - x_first) / x_step);
x_idx_max = round((upper_right_lon - x_first) / x_step);
y_idx_min = round((upper_right_lat - y_first) / y_step)+1;
y_idx_max = round((lower_left_lat - y_first) /  y_step)+1;

% Extract the region of interest for all time slices
region_data = timeseries_data(x_idx_min:x_idx_max, y_idx_min:y_idx_max, :);

% Calculate the sum of absolute displacements across all 39 time slices
sum_abs_displacement = sum(abs(region_data), 3);

% Mask the pixels where the sum of absolute displacements is greater than 20 cm
threshold_cm = 100;
mask = sum_abs_displacement * 100 > threshold_cm; % Convert to cm and apply threshold

% Initialize the average time series array
num_time_slices = size(timeseries_data, 3);
average_time_series = zeros(1, num_time_slices);

% Iterate through each time slice and calculate the average for non-masked values
for t = 1:num_time_slices
    % Extract the time slice for the region
    time_slice = region_data(:, :, t);
    
    % Apply the mask: set masked pixels to NaN
    time_slice(mask) = NaN;
    
    % Calculate the mean ignoring NaN values
    average_time_series(t) = mean(time_slice(:), 'omitnan');
end

% Convert average time series from meters to centimeters
average_time_series_cm = average_time_series * 100;


% Subtract the specific point time series from the average time series
resulting_time_series = point_time_series_cm - average_time_series_cm';


% Convert dates to MATLAB datetime objects
dates = h5read(filePath, '/date');
dates = cellstr(dates'); % Convert to cell array of strings
time_series_dates = datetime(dates, 'InputFormat', 'yyyyMMdd');

% Plot the resulting time series and average time series as subplots
figure;

% Upper Plot: Difference between specific point and average time series
subplot(2, 1, 1); % Create a 2-row, 1-column grid, and select the first plot
plot(time_series_dates, resulting_time_series, '-o', 'LineWidth', 1, 'Color', 'r','DisplayName', 'After Filter');
hold on; 
plot(time_series_dates, point_time_series_cm, '-o', 'LineWidth', 0.6, 'Color', 'b','DisplayName', 'Before Filter');
hold off; % Release the hold
xlabel('Date');
ylabel('Difference in Displacement (cm)');
title('Difference Between TS Displacement Before and After Applying Common Mode Filtering');
% grid on;
% Set the y-axis limits for both subplots
ylim([-3.5 3.5]);
legend('show', 'Location', 'northwest'); % Place legend in the upper right corner

% Lower Plot: Average displacement time series
subplot(2, 1, 2); % Select the second plot in the 2-row, 1-column grid
plot(time_series_dates, average_time_series_cm, '-o', 'LineWidth', 1, 'Color', 'k','DisplayName', 'Avg');
xlabel('Date');
ylabel('Average Displacement (cm)');
title('Average Time Series for Masked Region (No deformation area)');
grid on;
legend('show', 'Location', 'northwest'); % Place legend in the upper right corner

% % Link the y-axes of both plots to have the same y-axis scale
% linkaxes(findall(gcf, 'type', 'axes'), 'y');
% Save the figure in high quality
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperPosition', [0 0 35 25]); % Adjust this to stretch the figure horizontally
% set(gcf, 'PaperSize', [35 25]); % Match the paper size to the paper position for proper scaling
% print('TS_Des_p6_RP_ERA5.pdf', '-dpdf', '-r500'); % Save as PNG with 500 dpi resolution

% %------------------------ SAVE --------------------------------------------
% % TO SAVE in TXT file
% 
% % Convert the datetime array to the desired string format (YYYYMMDD)
% formatted_dates = datestr(time_series_dates, 'yyyymmdd');
% 
% % Combine dates and resulting time series into a cell array
% output_data = cell(length(formatted_dates), 2);
% output_data(:, 1) = cellstr(formatted_dates); % First column: formatted dates
% 
% output_data(:, 2) = num2cell(point_time_series_cm); % Second column: resulting time series
% % output_data(:, 2) = num2cell(resulting_time_series); % Second column: resulting time series
% 
% % Specify the output file name
% output_file = 'P5_mintpy_Desc.txt';
% 
% % Open the file for writing
% fileID = fopen(output_file, 'w');
% 
% % Write the header (optional)
% fprintf(fileID, 'Date\tResulting_Time_Series\n');
% 
% % Write the data to the file
% for i = 1:length(output_data)
%     fprintf(fileID, '%s\t%.4f\n', output_data{i, 1}, output_data{i, 2});
% end
% 
% % Close the file
% fclose(fileID);

% disp(['Data saved to ', output_file]);
%--------------------------------------------------------------------------
% Convert the datetime array to the desired string format (YYYYMMDD)
formatted_dates = datestr(time_series_dates, 'yyyymmdd');

% Combine dates and resulting time series into a cell array
output_data = cell(length(formatted_dates), 2);
output_data(:, 1) = cellstr(formatted_dates); % First column: formatted dates

% output_data(:, 2) = num2cell(point_time_series_cm); % Second column: resulting time series
output_data(:, 2) = num2cell(resulting_time_series); % Second column: resulting time series

% Specify the output file name
output_file = 'P5_mintpy_CMF_Desc.txt';

% Open the file for writing
fileID = fopen(output_file, 'w');

% Write the header (optional)
fprintf(fileID, 'Date\tResulting_Time_Series\n');

% Write the data to the file
for i = 1:length(output_data)
    fprintf(fileID, '%s\t%.4f\n', output_data{i, 1}, output_data{i, 2});
end

% Close the file
fclose(fileID);



%% 
% _________________________________________________________________________
% Histogram and map masked
%__________________________________________________________________________
% Get the number of slices
numSlices = size(region_data, 3);
% Loop through each slice and plot
for i = numSlices:numSlices
    % Create a new figure
    figure;

    % Extract the i-th slice
    data_slice = region_data(:, :, i) * 100; % Convert to cm

    % Plot the original data slice on the left
    subplot(1, 2, 1); % 1 row, 2 columns, 1st subplot
    imagesc(data_slice);
    colormap(jet);
    clim([-7 5]);
    colorbar;
    axis image; % Keep the aspect ratio equal
    title(['Original Data - Time Slice ', num2str(i)]);

    % Plot the masked region data on the right
    data_slice(mask) = NaN; % Apply mask to the data slice
    subplot(1, 2, 2); % 1 row, 2 columns, 2nd subplot
    imagesc(data_slice);
    colormap(jet);
    clim([-7 5]);
    colorbar;
    axis image;
    title(['Masked Region Data - Time Slice ', num2str(i)]);

    % Pause to allow viewing of the figure before proceeding to the next one
    pause(0.5); % Adjust this to your preference (0.5 seconds)
end


% Flatten the matrix to a vector
sum_abs_displacement_vector = sum_abs_displacement(:)*100;
% Plot the histogram
figure;
histogram(sum_abs_displacement_vector, 'BinWidth', 0.5); % Adjust BinWidth as needed
xlabel('Sum of Absolute Displacement (cm)');
ylabel('Frequency');
title('Histogram of Sum of Absolute Displacement for Selected Region');
grid on;

%%

%--------------------------------------------------------------------------
% To PLOT AVG for REF point of three h5 files 
% filePath = 'timeseries_D_R.h5';
% filePath = 'timeseries_ERA5_D_R.h5';
% filePath = 'Taftan_timeseries_Corrected.h5';



clc;
clear;
close all;

% Define file paths for all HDF5 files
filePaths = {'timeseries_D_R.h5', 'timeseries_ERA5_D_R.h5', 'timeseries_ERA5_demErr_D_R.h5'};

% Initialize cell arrays to store the results for each HDF5 file
average_time_series_cm_all = cell(1, numel(filePaths));
time_series_dates = [];

% Loop through each file and perform the analysis
for i = 1:numel(filePaths)
    % Load the HDF5 file
    filePath = filePaths{i};

    % Load the 'timeseries' dataset
    timeseries_data = h5read(filePath, '/timeseries');

    % Load metadata for coordinate calculations from the HDF5 file
    x_first = str2double(h5readatt(filePath, '/', 'X_FIRST'));
    x_step = str2double(h5readatt(filePath, '/', 'X_STEP'));
    y_first = str2double(h5readatt(filePath, '/', 'Y_FIRST'));
    y_step = str2double(h5readatt(filePath, '/', 'Y_STEP'));
    ncols = size(timeseries_data, 2); % Assuming data is [nrows x ncols x ntimes]
    nrows = size(timeseries_data, 1);

    % -------------------- Select POINT (e.g Ref. point)
    % Get the time series for the specific point  (REF Point)
    target_lat = 28.62662;
    target_lon = 61.14602;

    % Calculate indices for the target point
    target_x_idx = round((target_lon - x_first) / x_step);
    target_y_idx = round((target_lat - y_first) / y_step);

    % Define the bounding box (lower left and upper right corners)
    lower_left_lon = 61.11;
    lower_left_lat = 28.58;
    upper_right_lon = 61.15;
    upper_right_lat = 28.62;

    % Calculate indices for the bounding box
    x_idx_min = round((lower_left_lon - x_first) / x_step);
    x_idx_max = round((upper_right_lon - x_first) / x_step);
    y_idx_min = round((upper_right_lat - y_first) / y_step) + 1;
    y_idx_max = round((lower_left_lat - y_first) / y_step) + 1;

    % Extract the region of interest for all time slices
    region_data = timeseries_data(x_idx_min:x_idx_max, y_idx_min:y_idx_max, :);

    % Calculate the sum of absolute displacements across all time slices
    sum_abs_displacement = sum(abs(region_data), 3);

    % Mask the pixels where the sum of absolute displacements is greater than 20 cm
    threshold_cm = 180;
    mask = sum_abs_displacement * 100 > threshold_cm; % Convert to cm and apply threshold

    % Initialize the average time series array
    num_time_slices = size(timeseries_data, 3);
    average_time_series = zeros(1, num_time_slices);

    % Iterate through each time slice and calculate the average for non-masked values
    for t = 1:num_time_slices
        % Extract the time slice for the region
        time_slice = region_data(:, :, t);

        % Apply the mask: set masked pixels to NaN
        time_slice(mask) = NaN;

        % Calculate the mean ignoring NaN values
        average_time_series(t) = mean(time_slice(:), 'omitnan');
    end

    % Convert average time series from meters to centimeters
    average_time_series_cm = average_time_series * 100;

    % Store the results for plotting later
    average_time_series_cm_all{i} = average_time_series_cm;

    % Load dates (only once, as they are the same for all files)
    if isempty(time_series_dates)
        dates = h5read(filePath, '/date');
        dates = cellstr(dates'); % Convert to cell array of strings
        time_series_dates = datetime(dates, 'InputFormat', 'yyyyMMdd');
    end
end

% Find the global minimum and maximum y-axis values for consistent scaling
all_values = cell2mat(average_time_series_cm_all);
y_min = min(all_values(~isnan(all_values)));
y_max = max(all_values(~isnan(all_values)));
if y_min >= y_max
    y_min = y_min - 1;
    y_max = y_max + 1;
end

% Plot the average time series for all files in three subplots
figure;
for i = 1:numel(filePaths)
    subplot(3, 1, i); % Create a 3-row, 1-column grid, and select the appropriate plot
    plot(time_series_dates, average_time_series_cm_all{i}, '-o', 'LineWidth', 1, 'MarkerSize', 4, 'Color', 'k', 'DisplayName', 'Avg');
    xlabel('Date');
    ylabel('Average Displacement (cm)');
    title(['Average Time Series for Masked Region (No deformation area) (' filePaths{i} ')']);
    grid on;
    legend('show', 'Location', 'northwest'); % Place legend in the upper right corner
    ylim([y_min, y_max]); % Set consistent y-axis limits for all subplots
end

% Link the x-axes of all plots to have the same x-axis scale
linkaxes(findall(gcf, 'type', 'axes'), 'x');
% % Save the figure in high quality
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperPosition', [0 0 35 25]); % Adjust this to stretch the figure horizontally
% set(gcf, 'PaperSize', [35 25]); % Match the paper size to the paper position for proper scaling
% print('TS_avg_all_D_RP', '-dpng', '-r500'); % Save as PNG with 500 dpi resolution





%%
%--------------------------------------------------------------------------
%  PLOT TIME SERIES OF THREE h5 FILES 
%--------------------------------------------------------------------------

clc;
clear;
close all;

% Define file paths for all HDF5 files
filePaths = {'timeseries_D_R.h5', 'timeseries_ERA5_D_R.h5', 'timeseries_ERA5_demErr_D_R.h5'};

% Initialize cell arrays to store the results for each HDF5 file
resulting_time_series_all = cell(1, numel(filePaths));
point_time_series_cm_all = cell(1, numel(filePaths));
time_series_dates = [];

% Loop through each file and perform the analysis
for i = 1:numel(filePaths)
    % Load the HDF5 file
    filePath = filePaths{i};

    % Load the 'timeseries' dataset
    timeseries_data = h5read(filePath, '/timeseries');

    % Load metadata for coordinate calculations from the HDF5 file
    x_first = str2double(h5readatt(filePath, '/', 'X_FIRST'));
    x_step = str2double(h5readatt(filePath, '/', 'X_STEP'));
    y_first = str2double(h5readatt(filePath, '/', 'Y_FIRST'));
    y_step = str2double(h5readatt(filePath, '/', 'Y_STEP'));
    ncols = size(timeseries_data, 2); % Assuming data is [nrows x ncols x ntimes]
    nrows = size(timeseries_data, 1);

    % -------------------- Select POINT (e.g Ref. point)
    % Get the time series for the specific point  (REF Point)
    % target_lat = 28.601013;
    % target_lon = 61.133124;

    target_lat =  28.602198;
    target_lon =  61.135630;

    % Calculate indices for the target point
    target_x_idx = round((target_lon - x_first) / x_step);
    target_y_idx = round((target_lat - y_first) / y_step);

    % Extract the time series for the specific point
    point_time_series = timeseries_data(target_x_idx, target_y_idx, :);         % (Row(X) , Col(Y))
    point_time_series_cm = squeeze(point_time_series) * 100; % Convert to cm

    % Define the bounding box (lower left and upper right corners)
    lower_left_lon = 61.11;
    lower_left_lat = 28.58;
    upper_right_lon = 61.15;
    upper_right_lat = 28.62;

    % Calculate indices for the bounding box
    x_idx_min = round((lower_left_lon - x_first) / x_step);
    x_idx_max = round((upper_right_lon - x_first) / x_step);
    y_idx_min = round((upper_right_lat - y_first) / y_step) + 1;
    y_idx_max = round((lower_left_lat - y_first) / y_step) + 1;

    % Extract the region of interest for all time slices
    region_data = timeseries_data(x_idx_min:x_idx_max, y_idx_min:y_idx_max, :);

    % Calculate the sum of absolute displacements across all time slices
    sum_abs_displacement = sum(abs(region_data), 3);

    % Mask the pixels where the sum of absolute displacements is greater than 20 cm
    threshold_cm = 180;
    mask = sum_abs_displacement * 100 > threshold_cm; % Convert to cm and apply threshold

    % Initialize the average time series array
    num_time_slices = size(timeseries_data, 3);
    average_time_series = zeros(1, num_time_slices);

    % Iterate through each time slice and calculate the average for non-masked values
    for t = 1:num_time_slices
        % Extract the time slice for the region
        time_slice = region_data(:, :, t);

        % Apply the mask: set masked pixels to NaN
        time_slice(mask) = NaN;

        % Calculate the mean ignoring NaN values
        average_time_series(t) = mean(time_slice(:), 'omitnan');
    end

    % Convert average time series from meters to centimeters
    average_time_series_cm = average_time_series * 100;

    % Subtract the specific point time series from the average time series
    resulting_time_series = point_time_series_cm - average_time_series_cm';

    % Store the results for plotting later
    resulting_time_series_all{i} = resulting_time_series;
    point_time_series_cm_all{i} = point_time_series_cm;

    % Load dates (only once, as they are the same for all files)
    if isempty(time_series_dates)
        dates = h5read(filePath, '/date');
        dates = cellstr(dates'); % Convert to cell array of strings
        time_series_dates = datetime(dates, 'InputFormat', 'yyyyMMdd');
    end
end

% Find the global minimum and maximum y-axis values for consistent scaling
all_values = [cell2mat(resulting_time_series_all), cell2mat(point_time_series_cm_all)];
y_min = min(all_values(~isnan(all_values)));
y_max = max(all_values(~isnan(all_values)));
if y_min >= y_max
    y_min = y_min - 1;
    y_max = y_max + 1;
end
if y_min == y_max
    y_min = y_min - 1;
    y_max = y_max + 1;
end

% Plot the resulting time series for all files in three subplots
figure;
for i = 1:numel(filePaths)
    subplot(3, 1, i); % Create a 3-row, 1-column grid, and select the appropriate plot
    plot(time_series_dates, resulting_time_series_all{i}, '-o', 'LineWidth', 1, 'MarkerSize', 4, 'Color', 'b', 'DisplayName', 'After Filter');
    hold on;
    plot(time_series_dates, point_time_series_cm_all{i}, '-o', 'LineWidth', 1, 'MarkerSize', 4, 'Color', 'r', 'DisplayName', 'Before Filter');
    hold off;
    xlabel('Date');
    ylabel('Difference in Displacement (cm)');
    title(['Difference Between TS Displacement Before and After Applying Common Mode Filtering (' filePaths{i} ')']);
    grid on;
    legend('show', 'Location', 'northwest'); % Place legend in the upper right corner
    ylim([y_min, y_max]); % Set consistent y-axis limits for all subplots
end

% Link the x-axes of all plots to have the same x-axis scale
linkaxes(findall(gcf, 'type', 'axes'), 'x');
% % Save the figure in high quality
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperPosition', [0 0 35 25]); % Adjust this to stretch the figure horizontally
% set(gcf, 'PaperSize', [35 25]); % Match the paper size to the paper position for proper scaling
% print('TS_all_D_RP_p6', '-dpng', '-r500'); % Save as PNG with 500 dpi resolution



%##########################################################################



%%
%--------------------------------------------------------------------------
% TO SAVE in TXT file
%--------------------------------------------------------------------------


% Convert the datetime array to the desired string format (YYYYMMDD)
formatted_dates = datestr(time_series_dates, 'yyyymmdd');

% Combine dates and resulting time series into a cell array
output_data = cell(length(formatted_dates), 2);
output_data(:, 1) = cellstr(formatted_dates); % First column: formatted dates
output_data(:, 2) = num2cell(resulting_time_series); % Second column: resulting time series

% Specify the output file name
output_file = 'P0_summit_output_Des_CMF.txt';

% Open the file for writing
fileID = fopen(output_file, 'w');

% Write the header (optional)
fprintf(fileID, 'Date\tResulting_Time_Series\n');

% Write the data to the file
for i = 1:length(output_data)
    fprintf(fileID, '%s\t%.4f\n', output_data{i, 1}, output_data{i, 2});
end

% Close the file
fclose(fileID);

% disp(['Data saved to ', output_file]);
