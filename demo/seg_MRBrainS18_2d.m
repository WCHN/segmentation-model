% This script demonstrates using the trained segmentation-model to segment
% a subject from the MRBrainS18 segmentation challenge:
%     http://mrbrains18.isi.uu.nl/
%
% The output segmentation can be found in './output/segment/s/seg', along
% with other stuff.
% _________________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

clear; clc; close all;

% Add to the MATLAB path
addpath('..');

%--------------------------------------------------------------------------
% Data directory
%--------------------------------------------------------------------------

dir_data = './data-MRBrainS18-2d';

%--------------------------------------------------------------------------
% Set options (opt)
%--------------------------------------------------------------------------

% These two are mandatory (for now)
opt.dep.aux_toolbox  = '/path/to/auxiliary-functions';   % https://github.com/WTCN-computational-anatomy-group/auxiliary-functions
opt.dep.dist_toolbox = '/path/to/distributed-computing'; % https://github.com/WTCN-computational-anatomy-group/distributed-computing

opt.dir_output        = './output/';
opt.template.do       = false;
opt.verbose.level     = 1; % [0,1,2]
opt.seg.bg            = [7 8 9];
opt.template.bg_class = 7;

% Path to learnt model parameters
dir_model                 = './output/train/model';
opt.template.pth_template = fullfile(dir_model,'template.nii');
opt.gmm.pth_GaussPrior    = fullfile(dir_model,'GaussPrior.mat');
opt.gmm.pth_PropPrior     = fullfile(dir_model,'PropPrior.mat');

% For using multiple Gaussians per tissue
map          = containers.Map;
map('MRI')   = [1 2 3 4 5 5 6 6 7 7 8 8 8 9 9];
opt.dict.lkp = map; 

%--------------------------------------------------------------------------
% Segment subject(s)
%--------------------------------------------------------------------------

res = SegModel('segment',dir_data,opt);