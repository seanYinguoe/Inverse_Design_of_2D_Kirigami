function root = setup_project
%SETUP_PROJECT Add this paper's functions and examples, without saving paths.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root, 'src', 'kinematics')));
addpath(fullfile(root, 'examples'));
if ~isfolder(fullfile(root, 'results')), mkdir(fullfile(root, 'results')); end
end
