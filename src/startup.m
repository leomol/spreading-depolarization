% Download and unzip example dataset.
% Data will be downloaded once unless the downloads folder within the project folder is removed.

% 2023-09-01. Leonardo Molina
% 2024-09-25. Last modified.
setup();
function setup()
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'Chronux'));
    folder = fullfile(parent(root), 'downloads');
    target = fullfile(parent(root), 'download.zip');
    if exist(folder, 'dir') ~= 7
        fprintf('Downloading dataset to "%s"... ', target);
        source = 'https://github.com/leomol/spreading-depolarization-data/archive/refs/tags/v0.0.2.zip';
        filename = websave(target, source);
        fprintf('Done!\n');
        fprintf('Unzipping dataset"... ');
        unzip(filename, folder);
        fprintf('Done!\n');
    end
end

function path = parent(path)
    path = strsplit(path, filesep);
    path = fullfile(path{1:end - 1});
end