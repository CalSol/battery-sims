function setup()
%SETUP Add repository source, model, prototype, and diagnostic folders to the MATLAB path.

root = fileparts(mfilename('fullpath'));
folders = {'src', 'models', 'prototypes', 'tools'};

for k = 1:numel(folders)
    addpath(genpath(fullfile(root, folders{k})));
end

fprintf('battery-sims paths configured for this MATLAB session.\n');
end
