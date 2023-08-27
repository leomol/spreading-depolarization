% Usage:
%   First time ever:
%     -Install toolboxes: Curve Fitting, Signal Processing, Statistics and Machine Learning.
%     -Download and extract FPA from github.com/leomol/FPA, move it to Documents/MATLAB.
%   Every time MATLAB is reopened:
%     -run startup.m under Documents/MATLAB/FPA
%     -Run this script.
% 
% Position and FP data is assummed to start synchronously.
% Position data where DLC's prediction confidence is below a threshold is interpolated.
% Interpolation is not possible at the edges. These data is removed from both position and PF traces.

% 2022-12-07. Leonardo Molina.
% 2023-08-29. Last modified.

%% Configuration.
% Path to files containing fibre-photometry and dlc data (e.g. data/F46-doric.csv, data/F46-dlc.csv, ...)
fpFile = '../data/F21-doric.csv';
dlcFile = '../data/F21-dlc.csv';
[folder, session] = fileparts(fpFile);
% Number of pixels per cm.
ppcm = 15.40;

% Fiber photometry data.
frequency = 100;
baselineEpoch1 = [300, 600];
baselineEpoch2 = [4200, Inf];
thresholdEpochs = [200, 600];
sampleRange = [-Inf, Inf];
thresholdWindow = 5;
maxXCorrLag = 65;

% Columns corresponding to 465nm and 405nm.
data = loadData(fpFile);
time = data(:, 1);
lesion465 = data(:, 2);
lesion405 = data(:, 3);
nonlesion465 = data(:, 4);
nonlesion405 = data(:, 5);

% Position data from a DLC file.
% Arena size (cm).
arenaSize = [50, 50];
% Video frame rate.
framerate = 20;
% Minimum confidence level from DLC.
dlcThreshold = 0.90;
% Heat map saturation.
heatmapSaturation = 5;
heatmapBinCount = 20;

% General configuration.
configuration = struct();
configuration.resamplingFrequency = frequency;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 10.0;
configuration.peakSeparation = 0.5;
configuration.threshold = @(fpa, data) 4 * std(data - movmedian(data, round(thresholdWindow * frequency)));
configuration.thresholdEpochs = thresholdEpochs;
configuration.peakDetectionMode = 'prominence';
configuration.peakWindow = 2.5;
% Modify baseline correction so that bleaching estimated from reference is subtracted from both signal and reference.
configuration.signalCorrection = @(fpa, data) data - fpa.referenceBaseline;
configuration.referenceCorrection = @(fpa, data) data - fpa.referenceBaseline;
% Use a portion of the data to compute modified z-score.
configuration.f0 = {@median, baselineEpoch1};
configuration.f1 = {@mad, baselineEpoch1};
% Define epochs for analyses (FP & behaviour); Ipsi = IL stroke event; CL or contra = CL seizure/spreading depolarization
configuration.conditionEpochs = {'Baseline', baselineEpoch1, 'PreIpsi', [501, 801], 'Ipsi', [751, 851], 'PostIpsi', [801, 1101], 'PreCL1', [760, 830], 'CL', [1385, 1485], 'PostCL1', [2330, 2480], 'Final', [2100,2400]};
% Use a portion of the data to model and correct for photo-bleaching. Signal uses only baselineEpoch1 while reference uses baselineEpoch1 and baselineEpoch2.
configuration.baselineEpochs = {baselineEpoch1, [baselineEpoch1, baselineEpoch2]};

% Load behavioral data.
if isfile(dlcFile)
    dlc = loadDLC(dlcFile);
    
    % Choose x and y for body position.
    x = dlc.bodypart2_x / ppcm;
    y = dlc.bodypart2_y / ppcm;
    p = dlc.bodypart2_p;
    nFrames = size(dlc, 1);
    
    % Filter.
    cutoff = framerate / 2;
    order = 9;
    db = 20;
    [b, a] = cheby2(order, db, cutoff / framerate / 2, 'low');
    k2 = p > dlcThreshold;
    x2 = filtfilt(b, a, x(k2));
    y2 = filtfilt(b, a, y(k2));
    t1 = (0:nFrames - 1)' / framerate;
    t1 = t1(k2);
    
    % Clip traces traces. Assume that behavior starts synchronously with FP.
    minTime = max(t1(1), time(1));
    maxTime = min(t1(end), time(end));
    maxTime = floor(maxTime * frequency) / frequency;
    t2 = (minTime:1 / frequency:maxTime)';
    x2 = interp1(t1, x2, t2);
    y2 = interp1(t1, y2, t2);
    
    k1 = time >= minTime & time <= maxTime;
    time = time(k1);
    lesion465 = lesion465(k1);
    lesion405 = lesion405(k1);
    nonlesion465 = nonlesion465(k1);
    nonlesion405 = nonlesion405(k1);
else
    fprintf(2, '"%s" not found.\n', dlcFile);
    t2 = time;
end

% Call FPA with given configuration.
lesionFpa = FPA(time, lesion465, lesion405, configuration);
cellfun(@warning, lesionFpa.warnings);
nonLesionFpa = FPA(time, nonlesion465, nonlesion405, configuration);
cellfun(@warning, nonLesionFpa.warnings);

% Plot all fpa related figures.
lesionFpa.plotTrace();
lesionFpa.plotStatistics();
nonLesionFpa.plotTrace();
nonLesionFpa.plotStatistics();

%% Option 1 - Save all FPA data.
lesionFpa.export([session, ' - lesion']);
nonLesionFpa.export([session, ' - non-lesion']);

%% Option 2 - Save dff with lower sampling rate for plotting.
samplingRate = 10;
filename = fullfile(folder, sprintf('%s - dff plot - lesion.csv', session));
export(filename, lesionFpa.time, lesionFpa.dff, samplingRate);

%% Locomotion trace.
figure();
ids1 = time2id(t1(1:end-1), sampleRange);
ids2 = time2id(t2(1:end-1), sampleRange);
ax = NaN(1, 2);
ax(1) = subplot(1, 2, 1);
plot(x(ids1), y(ids1));
title('Raw');
xlabel('x (cm)');
ylabel('y (cm)');
axis('square');
ax(2) = subplot(1, 2, 2);
plot(x2(ids2), y2(ids2));
title('Processed');
xlabel('x (cm)');
axis('square');
linkaxes(ax);
axis(ax(1), 'tight');

%% Locomotion over time.
delta = [0; sqrt(diff(x2) .^ 2 + diff(y2) .^ 2)];
figure();
plot(t2, cumsum(delta));
axis('tight');
xlabel('time (s)');
ylabel('Locomotion (cm)');

% Total locomotion for individual epochs.
nEpochs = numel(configuration.conditionEpochs) / 2;
labels = configuration.conditionEpochs(1:2:end);
locomotion = zeros(nEpochs, 1);
for e = 1:nEpochs
    range = configuration.conditionEpochs{2 * e - 0};
    k = time2id(t2, range);
    locomotion(e) = sum(delta(k));
end
figure()
bar(locomotion)
xticklabels(labels);
xlabel('Epoch');
ylabel('Locomotion (cm)');
xtickangle(45);
grid('on');
grid('minor');

%% Locomotion heatmap.
ids = time2id(t2(1:end-1), sampleRange);
xlims = [0, arenaSize(1)] + min(x);
ylims = [0, arenaSize(2)] + min(y);
xEdges = linspace(xlims(1), xlims(2), heatmapBinCount);
yEdges = linspace(ylims(1), ylims(2), heatmapBinCount);
[xGrid, yGrid] = meshgrid(xEdges, yEdges);
xGrid = xGrid(1:end - 1, 1:end - 1);
yGrid = yGrid(1:end - 1, 1:end - 1);
intensity = histcounts2(y2(ids), x2(ids), yEdges, xEdges);
figure();
pcolor(xGrid, yGrid, intensity);
shading('interp');
clim([0, max(prctile(intensity(:), 100 - heatmapSaturation), 1)]);
axis('square');
xlabel('x (cm)');
ylabel('y (cm)');

%% Speed.
speed = sqrt(x2 .^ 2 + y2 .^ 2);
speed = diff(speed) ./ diff(t2);
speed = [speed; speed(end)];
filteredSpeed = movmedian(speed, 2 * frequency);
figure();
plot(t2, abs(filteredSpeed));
xlabel('time (s)');
ylabel('v (cm/s)');
axis('tight');

%% Cross-correlations.
ids = time2id(t2(1:end-1), sampleRange);
[xc1,    ~] = xcorr(nonLesionFpa.dff(ids), speed(1:numel(nonLesionFpa.dff(ids))), round(maxXCorrLag * frequency), 'normalized');
[xc2, lags] = xcorr(lesionFpa.dff(ids), speed(1:numel(lesionFpa.dff(ids))), round(maxXCorrLag * frequency), 'normalized');
figure();
subplot(2, 2, 1);
plot(lags / frequency, xc1, 'b');
axis('tight');
subplot(2, 2, 3);
plot(lags / frequency, xc2, 'r');
axis('tight');
subplot(2, 2, [2, 4]);
hold('all');
plot(lags / frequency, xc1, 'b', 'DisplayName', 'XCorr of speed to non-lesion signal');
plot(lags / frequency, xc2, 'r', 'DisplayName', 'XCorr of speed to lesion signal');
mn = min([xc1; xc2]);
mx = max([xc1; xc2]);
plot([0, 0], [mn, mx], 'Color', [0.50, 0.50, 0.50], 'LineStyle', '--', 'Marker', 'none', 'HandleVisibility', 'off');
axis('tight');
legend('show');

function export(output, time, x, frequency, header)
    [p, q] = rat(frequency * median(diff(time)));
    [x, time] = resample(x, time, frequency, p, q);
    fid = fopen(output, 'w');
    if nargin == 5
        fprintf(fid, [header, '\n']);
    end
    fprintf(fid, '%f, %f\n', [time, x]');
    fclose(fid);
end