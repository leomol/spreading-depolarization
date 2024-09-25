% 2023-11-30. LM added spectrograms using the multi-tapers method.
% 2024-09-16. AB added spectrograms using cosine-bell windows.
% 2024-09-25. Last modified.

%% User configuration.

% Define input and output data.
data_folder = '../downloads/spreading-depolarization-data-0.0.2/data/';
fp_file = fullfile(data_folder, 'ExampleStrokeMK801-GCaMP.csv');
lfp_file = fullfile(data_folder, 'ExampleStrokeMK801-LFP.csv');
output_filename = 'outputs/ExampleStrokeMK801';

% Define time range for power analysis of seizure in IL vs CL LFP.
seizure_period = [1575, 1635];

% Define the time windows (2) for each segment (duration in seconds) for CL_LFP PSD calculation (FFT, cosine-bell).
time_windows = [300, 360; 1575, 1635];
segment_duration = 60; % Duration of each segment in seconds.
window_size = 256;
overlap = 0.9375; % 93.75% overlap.

% Parameters for time windows of interest (BL vs post-photothrombosis) for suprathreshold peak detection.
baseline_period = [300, 360];  % Baseline periods in seconds.
analysis_period = [601, 2400]; % Analysis period in seconds.

% Threshold factor for amplitude-based suprathreshold peak detection:
threshold_factor = 4;

% Define settings for spectrogram calculation.
mts_time_bandwidth_product = 20;
mts_n_tapers = 9;
mts_window_size = 10;
mts_window_step = 1;

% Whether to save data files or not.
save_data = true;

%% Load data.

% Add FPA and Chronux paths.
try
    run('FPA-2.0.0/src/startup.m');
catch
end

% Spectrogram plots require Chronux library.
if exist('mtspecgramc', 'file') ~= 2
    error('Chronux library required!');
end

% Load FP data (from FPA analysis output) with columns for time, IL, CL.
data = CSV.load(fp_file);
time_G = data(:, 1);
IL_G = data(:, 2);
CL_G = data(:, 3);

% Load LFP data (from .csv converted from WINDAQ file) with columns for time, IL, CL.
data = CSV.load(lfp_file);
time_LFP = data(:, 1);
IL_LFP = data(:, 3); % IL LFP column.
CL_LFP = data(:, 2); % CL LFP column.

% Create output folder.
folder = fileparts(output_filename);
d = warning('QUERY', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(folder);
warning(d.state, 'MATLAB:MKDIR:DirectoryExists');

%% Run spectrogram plot for above signal.
sampling_frequency = 1 / mean(diff(time_LFP)); % 256.14235;

figure('Name', 'Unfiltered CL Spectrogram');
subplot(4, 1, 1);
plot(time_G, IL_G, 'k');
ylabel('IL GCaMP6f (DF/F)');
title('Ipsilesional GCaMP6f');

subplot(4, 1, 2);
plot(time_G, CL_G, 'r');
ylabel('CL GCaMP (DF/F)');
title('Contralesional GCaMP6f');

subplot(4, 1, 3);
plot(time_LFP, CL_LFP, 'r');
ylabel('CL LFP (µV)')
title('Contralesional LFP');

subplot(4,1,4);
params = struct('Fs', sampling_frequency, 'tapers', [mts_time_bandwidth_product, mts_n_tapers]);
[psd_estimates, tics, frequencies] = mtspecgramc(CL_LFP, [mts_window_size, mts_window_step], params);
imagesc(tics, frequencies, 10 * log10(psd_estimates)');
axis('xy');
colormap('jet');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Contralesional LFP Spectrogram');

% Link x-axis.
axs = findobj(gcf(), 'type', 'axes');
linkaxes(axs, 'x');
xlim(axs(1), [time_G(1), time_G(end)]);
set(axs(2:end), 'XTick', []);

%% Plot LFP spectrogram for both sides.
figure('Name', 'Bilateral LFP Spectrogram');
subplot(2, 1, 1);
params = struct('Fs', sampling_frequency, 'tapers', [mts_time_bandwidth_product, mts_n_tapers]);
[psd_estimates, tics, frequencies] = mtspecgramc(IL_LFP, [mts_window_size, mts_window_step], params);
imagesc(tics, frequencies, 10 * log10(psd_estimates)');
axis('xy');
colormap('jet');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Ipsilesional LFP Spectrogram');

subplot(2, 1, 2);
params = struct('Fs', sampling_frequency, 'tapers', [mts_time_bandwidth_product, mts_n_tapers]);
[psd_estimates, tics, frequencies] = mtspecgramc(CL_LFP, [mts_window_size, mts_window_step], params);
imagesc(tics, frequencies, 10 * log10(psd_estimates)');
axis('xy');
colormap('jet');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Contralesional LFP Spectrogram');

% Link x-axis.
axs = findobj(gcf(), 'type', 'axes');
linkaxes(axs, 'x');
set(axs(2:end), 'XTick', []);
 
%% Identify and export PSD using FFT with cosine-bell for specified times.

lfp_PSD = CL_LFP; % Select whether you want IL or CL LFP analysed for PSD calculation.

overlap_samples = floor(window_size * overlap);

% Convert time windows to sample indices.
segment_samples = segment_duration * sampling_frequency;
time_windows_samples = round(time_windows * sampling_frequency);

% Preallocate cell array to store PSD data.
psd_data = cell(size(time_windows, 1), 1);

% Loop through each time window and compute PSD using FFT with a cosine-bell window.
for idx = 1:size(time_windows, 1)
    start_sample = time_windows_samples(idx, 1) + 1;
    end_sample = time_windows_samples(idx, 2);
    segment = lfp_PSD(start_sample:end_sample);
    
    % Apply cosine-bell window.
    cosine_bell_window = 1 - cos(2 * pi * (0:(window_size - 1))' / (window_size - 1));
    
    % Compute FFT with overlap.
    [S, F] = spectrogram(segment, cosine_bell_window, overlap_samples, window_size, sampling_frequency, 'yaxis');
    
    % Compute PSD.
    Pxx = abs(S) .^ 2 / (sampling_frequency * sum(cosine_bell_window .^ 2));
    % Convert to dB and average over segments.
    psd_data{idx} =  10 * log10(mean(Pxx, 2));
end

% Plot the PSD for each time window.
figure();
hold('on');
for idx = 1:length(psd_data)
    plot(F, psd_data{idx}, 'DisplayName', sprintf('Time Window %d-%d s', time_windows(idx, 1), time_windows(idx, 2)));
end
hold('off');
xlim([0, 80]); % Frequency range from 0 to 80 Hz.
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density');
legend('show');

% Generate table containing PSD of each time window.
if save_data == true
    % Combine frequency values with PSD data.
    psd_table_data = cell(length(F), 3);
    psd_table_data(:, 1) = num2cell(F);
    
    for idx = 1:length(psd_data)
        psd_table_data(:, idx + 1) = num2cell(psd_data{idx});
    end
    
    % Create a table with PSD data.
    psd_table = cell2table(psd_table_data, 'VariableNames', {'Frequency', sprintf('baseline %d-%d s', time_windows(1, 1), time_windows(1, 2)), sprintf('seizure %d-%d s', time_windows(2, 1), time_windows(2, 2))});
    
    % Export table to CSV file.
    csv_file = sprintf('%s_CLpsd_data.csv', output_filename); % Specify your desired file name based on session.
    writetable(psd_table, csv_file);
end

%% Low-pass filter parameters (optional) prior to peak detection.

use_filter = true; % Set to true to apply low-pass filter.
cutoff_frequency = 20; % Cutoff frequency in Hz for the low-pass filter.

% Design low-pass filter if filtering is enabled.
if use_filter
    % Design a low-pass Butterworth filter
    [b, a] = butter(4, cutoff_frequency / (sampling_frequency / 2), 'low');
    % Apply the filter to both LFP signals
    IL_LFP = filtfilt(b, a, IL_LFP);
    CL_LFP = filtfilt(b, a, CL_LFP);
end

%% Identify mean and SD of baseline LFP.
% Calculate indices for baseline and analysis periods.
baseline_indices = round(baseline_period(1) * sampling_frequency) + 1 : round(baseline_period(2) * sampling_frequency);
analysis_indices = round(analysis_period(1) * sampling_frequency) + 1 : round(analysis_period(2) * sampling_frequency);

% Calculate baseline mean and standard deviation for IL_LFP.
baselineIL = IL_LFP(baseline_indices);
mean_baselineIL = mean(baselineIL);
std_baselineIL = std(baselineIL);

% Calculate baseline mean and standard deviation for CL_LFP.
baselineCL = CL_LFP(baseline_indices);
mean_baselineCL = mean(baselineCL);
std_baselineCL = std(baselineCL);

%% Label all timepoints where LFP > Threshold SD value.

% Identify timepoints in analysis period where lfp1 > threshFactor * std_baseline1.
IL_LFP_analysis = IL_LFP(analysis_indices);
thresholdIL = mean_baselineIL + threshold_factor * std_baselineIL;
exceeding_indicesIL = find(IL_LFP_analysis > thresholdIL | IL_LFP_analysis < -thresholdIL);
exceeding_timesIL = (exceeding_indicesIL - 1) / sampling_frequency + analysis_period(1);

% Identify timepoints in analysis period where lfp2 > threshFactor * std_baseline2.
CL_LFP_analysis = CL_LFP(analysis_indices);
thresholdCL = mean_baselineCL + threshold_factor * std_baselineCL;
exceeding_indicesCL = find(CL_LFP_analysis > thresholdCL | CL_LFP_analysis < -thresholdCL);
exceeding_timesCL = (exceeding_indicesCL - 1) / sampling_frequency + analysis_period(1);


%% Identify ONLY PEAKS where LFP > or < threshold factor (peaks defined by ThreshFactor).

% Identify positive peaks in analysis period (IL).
IL_LFP_analysis = IL_LFP(analysis_indices);
thresholdIL_pos = mean_baselineIL + threshold_factor * std_baselineIL;
[peak_valsIL_pos, peak_indicesIL_pos] = findpeaks(IL_LFP_analysis, 'MinPeakHeight', thresholdIL_pos);
peak_timesIL_pos = (peak_indicesIL_pos - 1) / sampling_frequency + analysis_period(1);

% Identify negative peaks in analysis period where (IL).
thresholdIL_neg = mean_baselineIL - threshold_factor * std_baselineIL;
[peak_valsIL_neg, peak_indicesIL_neg] = findpeaks(-IL_LFP_analysis, 'MinPeakHeight', -thresholdIL_neg);
peak_valsIL_neg = -peak_valsIL_neg; % Convert back to negative peaks
peak_timesIL_neg = (peak_indicesIL_neg - 1) / sampling_frequency + analysis_period(1);

% Identify positive peaks in analysis period (CL).
CL_LFP_analysis = CL_LFP(analysis_indices);
thresholdCL_pos = mean_baselineCL + threshold_factor * std_baselineCL;
[peak_valsCL_pos, peak_indicesCL_pos] = findpeaks(CL_LFP_analysis, 'MinPeakHeight', thresholdCL_pos);
peak_timesCL_pos = (peak_indicesCL_pos - 1) / sampling_frequency + analysis_period(1);

% Identify negative peaks in analysis period (CL).
thresholdCL_neg = mean_baselineCL - threshold_factor * std_baselineCL;
[peak_valsCL_neg, peak_indicesCL_neg] = findpeaks(-CL_LFP_analysis, 'MinPeakHeight', -thresholdCL_neg);
peak_valsCL_neg = -peak_valsCL_neg; % Convert back to negative peaks
peak_timesCL_neg = (peak_indicesCL_neg - 1) / sampling_frequency + analysis_period(1);

% Plot GCaMP6f (ipsi & contra).
figure();
subplot(3, 1, 1);
plot(time_G, IL_G, 'k', time_G, CL_G, 'r');
xlabel('Time (s)');
ylabel('∆F/F');
title('GCaMP');
legend('ipsi', 'contra');

% Plot LFP1 (ipsi).
subplot(3, 1, 2);
plot(time_LFP, IL_LFP, 'k');
hold('on');
plot(peak_timesIL_pos, peak_valsIL_pos, 'go', peak_timesIL_neg, peak_valsIL_neg, 'bo');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
title('Ipsi LFP with Exceeding Peaks');
legend('Ipsi', 'Exceeding Peaks');

% Plot LFP2 (contra).
subplot(3, 1, 3);
plot(time_LFP, CL_LFP, 'r');
hold('on');
plot(peak_timesCL_pos, peak_valsCL_pos, 'go', peak_timesCL_neg, peak_valsCL_neg, 'bo');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
title('Contra LFP with Exceeding Peaks');
legend('Contra', 'Exceeding Peaks');

% Link x-axis for all subplots.
linkaxes(findall(gcf, 'type', 'axes'), 'x');

%% Plot power over time for LFP signals.

% Define frequency bands.
freq_bands = [0, 4; 4, 8; 8, 12; 12, 30; 30, 80; 0, 80];  % Delta, Theta, Alpha, Beta, Gamma, All.
band_names = {'Delta (0-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-12 Hz)', 'Beta (12-30 Hz)', 'Gamma (30-80 Hz)', 'All (0-80 Hz)'};

% Initialize power matrix.
num_bands = size(freq_bands, 1);
power_IL = cell(num_bands, 1);
power_CL = cell(num_bands, 1);

% Parameters for STFT.
window_size = 256;
noverlap = window_size / 2;
nfft = 512;

for band = 1:num_bands
    % Compute the STFT for the current frequency band.
    [~, f, t_stft, P1] = spectrogram(IL_LFP, window_size, noverlap, nfft, sampling_frequency);
    [~, ~, ~, P2] = spectrogram(CL_LFP, window_size, noverlap, nfft, sampling_frequency);
    
    % Find indices corresponding to the current frequency band.
    freq_idx = find(f >= freq_bands(band, 1) & f <= freq_bands(band, 2));
    
    % Compute the average power in the current frequency band over time.
    power_IL{band} = 10 * log10(mean(P1(freq_idx, :), 1));
    power_CL{band} = 10 * log10(mean(P2(freq_idx, :), 1));
end


%% Plot the LFP & GCaMP with average power in common frequency bands (as well as overall power of signal).

% Plot GCaMP (IL and CL).
figure();
subplot(num_bands + 2, 1, 1); 
plot(time_G, IL_G, 'k', time_G, CL_G, 'r');
xlabel('Time (s)');
ylabel('∆F/F');
title('GCaMP');
legend('ipsi', 'contra');

% Plot LFP (IL and CL).
subplot(num_bands + 2, 1, 2);
plot(time_LFP, IL_LFP, 'k', time_LFP, CL_LFP, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('LFP Signals');
legend('ipsi', 'contra');

% Plot the average power over time for each frequency band.
for band = 1:num_bands
    nTstft = t_stft; 
    subplot(num_bands + 2, 1, band + 2);
    plot(nTstft, power_IL{band}, 'k', nTstft, power_CL{band}, 'r');
    xlabel('Time (s)');
    ylabel('Power (dB)');
    title(['Average Power - ' num2str(freq_bands(band, 1)) '-' num2str(freq_bands(band, 2)) ' Hz']);
    legend('ipsi', 'contra');
end

% Link x-axis for all subplots.
linkaxes(findall(gcf, 'type', 'axes'),'x');

%% Write table for average power in frequency bands over time (contra).

% Create an array to hold power data for all bands.
all_power_data2 = zeros(length(nTstft), num_bands);

% Populate the array with power data from each band.
for band = 1:num_bands
    all_power_data2(:, band) = power_CL{band};
end

% Create a table with time and power data.
cl_freq_tab = array2table(all_power_data2, 'VariableNames', band_names);
cl_freq_tab.Time = nTstft'; % Add the time vector as the first column.
cl_freq_tab = movevars(cl_freq_tab, 'Time', 'Before', 1); % Move 'Time' to be the first column.

if save_data == true
    writetable(cl_freq_tab, sprintf('%s-clFreqTab.csv', output_filename));
end

%% Create table of contra STFT power with normalization to max power at specific window around seizure.

% Find the indices corresponding to the desired time range.
time_indices = (cl_freq_tab.Time >= seizure_period(1)) & (cl_freq_tab.Time <= seizure_period(2));

% Filter the table to only include the specified time range.
filtered_cl_freq_tab = cl_freq_tab(time_indices, :);

% Normalize the power values for each band within the filtered range.
for band = 1:num_bands
    % Determine column name for the current band.
    column_name = band_names{band};

    % Calculate maximum power within the filtered time range.
    max_power2 = max(filtered_cl_freq_tab.(column_name)); 

    % Create new column name with "norm_" prefix.
    norm_column_name = ['norm_' column_name];

    % Normalize each power value to max CL power within the time range if max_power is not zero.
    if max_power2 == 0
        % If max_power is zero, set normalized values to zero.
        filtered_cl_freq_tab.(norm_column_name) = zeros(size(filtered_cl_freq_tab.(column_name)));
        disp('Max Power = 0 during time window');
    else
        filtered_cl_freq_tab.(norm_column_name) = filtered_cl_freq_tab.(column_name) / max_power2;
    end
end

if save_data == true
    writetable(filtered_cl_freq_tab, sprintf('%s-clWindowFreq.csv', output_filename));
end

%% Write table for average power in frequency bands over time (ipsi).

% Create an array to hold power data for all bands.
all_power_data = zeros(length(nTstft), num_bands);

% Populate the array with power data from each band.
for band = 1:num_bands
    all_power_data(:, band) = power_IL{band};
end

% Create a table with time and power data.
il_freq_tab = array2table(all_power_data, 'VariableNames', band_names);
il_freq_tab.Time = nTstft'; % Add the time vector as the first column.
il_freq_tab = movevars(il_freq_tab, 'Time', 'Before', 1); % Move 'Time' to be the first column.

if save_data == true
    writetable(il_freq_tab, sprintf('%s-ilFreqTab.csv', output_filename));
end

%% Create table of ipsi STFT power with normalization to max power at specific window around EE.

% Filter the table to only include the specified time range.
filtered_il_freq_tab = il_freq_tab(time_indices, :);

% Normalize the power values for each band within the filtered range.
for band = 1:num_bands
    % Determine column name for the current band.
    column_name = band_names{band};
    
    % Calculate maximum power within the filtered time range.
    max_power1 = max(filtered_il_freq_tab.(column_name)); 
    
    % Create new column name with "norm_" prefix
    norm_column_name = ['norm_' column_name];
    
    % Normalize each power value within the time range to Max CL power if max_power is not zero.
    if max_power2 == 0
        % If max_power is zero, set normalized values to zero.
        filtered_il_freq_tab.(norm_column_name) = zeros(size(filtered_il_freq_tab.(column_name)));
    elseif max_power1 <= max_power2
        filtered_il_freq_tab.(norm_column_name) = filtered_il_freq_tab.(column_name) / max_power2;
    elseif max_power1 > max_power2
        filtered_il_freq_tab.(norm_column_name) = filtered_il_freq_tab.(column_name) / max_power1;
    end
end

if save_data == true
    writetable(filtered_il_freq_tab, sprintf('%s-ilWindowFreq.csv', output_filename));
end