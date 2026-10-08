function root = setup_project
%SETUP_PROJECT Add only this project's source and configuration folders.
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'config'),genpath(fullfile(root,'src')));
if ~isfolder(fullfile(root,'results')), mkdir(fullfile(root,'results')); end
end
