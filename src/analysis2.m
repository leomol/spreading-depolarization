% 2023-08-29. Wilten Nicola.
% 2023-08-29. Last modified.

%% Setup.
cons = {'sitting', 'rearing', 'moving', 'grooming'};
files = dir('../data/F*.xlsx');

% Load xlsx files in a specified folder.
nFiles = numel(files);
filenames = {files.name};
paths = fullfile({files.folder}, filenames);
[time, ipsi, contr, state] = data_loader(paths, false, cons);

%% Process all animals to locate peak derivative.  
dt = time(2) - time(1);  % time step/inverse of sample frequency. 
xx = 10 / dt; % smoothing window. 
Q0 = 0; % initialize the behavioural average.
S0 = 0; % initialize the signal average.
t0 = -300 + dt:dt:300; % center data  around 0 from -300 to 300.  0 corresponds to the peak derivative. 
store_der = NaN(1, nFiles);
ipsi_time = NaN(1, nFiles);
colormap('hot');
for k = 1:nFiles
    yq = smooth(ipsi(:, k), xx, 'lowess');
    der = (yq(2:end) - yq(1:end - 1)) ./ (time(2:end) - time(1:end-1));
    der(time(1:end - 1) < 5) = 0; % pad the derivative at the edges to avoid boundary effects. 
    der(time(1:end - 1) > 2400) = 0;
    [mx,ix] = max(der);  % find the  maximum. 
    
    store_der(k) = mx; 
    ipsi_time(k) = time(ix);
    
    time_shift = time - ix * dt; % shift the time window.
    i0 = find(time_shift > -300 & time_shift < 300); 
    Q0 = Q0 + squeeze(state(i0, :, k));
    S0 = S0 + ipsi(i0, k);
    
    % Location of peak derivatives.
    figure(1);
    plot(time, yq / max(yq) + k);
    hold('on');
    plot(time(ix), yq(ix) / max(yq) + k, 'r.');
    title('Location of Derivative Peak');
    xlabel('Time');
    ylabel('Recording Index');
    xlim([0, max(time)]);
    
    % Peak derivatives.
    figure(2);
    subplot(4, 1, 1);
    plot(time(1:end-1) - ix * dt, der);
    hold('on');
    xlim([-300, 300]);
    title('Max Derivative');
    subplot(4, 1, 2);
    plot(time - ix * dt, ipsi(:, k));
    hold('on');
    xlim([-300, 300]);
    title('Ipsilateral (Max D Aligned)');
    
    % Behavior and peak derivatives.
    figure(3);
    subplot(nFiles, 2, 2 * k - 1)
    imagesc(-300-dt:dt:300, 1:4, squeeze(state(i0, :, k))');
    hold('on');
    plot(time(1:end-1) - ix * dt, 4 * der / max(der), 'r');
    plot(time(1:end-1) - ix * dt, 4 * contr(1:end - 1, k) / max(contr(:, k)), 'g', 'LineWidth', 2);
    xlim([-300, 300]);
    yticks([1, 2, 3, 4]);
    yticklabels(cons);
    set(gca, 'ydir', 'normal');
    colormap('hot');
    ylim([0, 4.5]);
    filename = filenames{k};
    title(filename(1:end-5));
    subplot(nFiles, 2, 2 * k);
    plot(time - ix * dt, ipsi(:, k) / max(ipsi(:, k)));
    hold('on');
    plot(time(1:end - 1) - ix * dt, der / max(der) + 0.5, 'r');
    xlim([-300, 300]);
    if k == nFiles
        legend('Ipsi', 'Derivative');
    end
end

% Normalize the means.
Q0 = Q0 / nFiles; 
S0 = S0 / nFiles;

% Plot the means.
figure(2);
subplot(4, 1, 3);
plot(-300 + dt:dt:300, Q0);
title('Behavioural Average');
legend(cons);
subplot(4, 1, 4);
plot(-300 + dt:dt:300, S0);
title('Ipsi Average');

function [time, ipsi, contr, state, dur] = data_loader(paths, anest, cons)
    anest = nargin < 2 | anest;
    if nargin < 3
        cons = {'sitting', 'rearing', 'moving', 'grooming'};
    end
    nCons = numel(cons);
    nFiles = numel(paths); 
    
    for j = 1:nFiles 
        num1 = xlsread(paths{j}, 1);
        time = num1(:, 1);
        ntime = length(time); 
        if j == 1 
        state = zeros(ntime, nCons, nFiles); 
        dur = [];
        ipsi = zeros(ntime, nFiles);
        contr = zeros(ntime, nFiles);
        end
        
        ipsi(1:length(num1(:, 2)), j) = num1(:, 2);
        if size(num1,2) > 2 
            contr(:,j) = num1(:, 3);
        else 
            contr(1:length(time), j) = zeros(size(time));
        end
        dt = time(2) - time(1);
        
        if ~anest
            [num2, txt2] = xlsread(paths{j}, 2);
            m = 0; 
            for k = 1:nCons
                i0 = find(strcmp(txt2(:, 1), cons{k})==1);
                num2(1) = dt; 
                start_ = num2(i0 - 1);
                end_ = num2(i0) - dt;
                for z = 1:length(end_)
                    m = m + 1; 
                    state(round(start_(z) / dt):round(end_(z) / dt), k, j) = 1; 
                    dur(j, m, :) = [start_(z), end_(z) - start_(z), k];
                end
            end
        end
    end
end