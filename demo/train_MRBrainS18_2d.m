% This script demonstrates using the segmentation-model code to learn from
% MR scans part of the MRBrainS18 segmentation challenge:
%     http://mrbrains18.isi.uu.nl/
%
% The code learns a template, hyper-parameters of a Gaussian mxiture model
% for intensities, and hyper-parameters of a dirichlet prior on the
% proportions of tissues. All these are stored in './output/train/model'.
%
% _________________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

clear; clc;

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
% opt.dep.aux_toolbox  = '/path/to/auxiliary-functions';   % https://github.com/WTCN-computational-anatomy-group/auxiliary-functions
% opt.dep.dist_toolbox = '/path/to/distributed-computing'; % https://github.com/WTCN-computational-anatomy-group/distributed-computing
opt.dep.aux_toolbox  = '/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions';   % https://github.com/WTCN-computational-anatomy-group/auxiliary-functions
opt.dep.dist_toolbox = '/home/mbrud/dev/mbrud/code/matlab/distributed-computing'; % https://github.com/WTCN-computational-anatomy-group/distributed-computing

% Directory for learnt model (and temp files)
opt.dir_output = './output';

% Template options
opt.template.do = true;
opt.template.K  = 10;

% % For mapping labels to tissue classes
map                = containers.Map;
map('MRBrainS18')  = {6,7,8,9,5,4,[]}; % 1.CGM 2.BG 3.WM 4.WMH 5.CSF 6.VEN
opt.gmm.labels.cm  = map;
opt.gmm.labels.use = true;

% For using multiple Gaussians per tissue
map          = containers.Map;
map('MRI')   = repelem(1:opt.template.K,2); 
opt.dict.lkp = map;

%--------------------------------------------------------------------------
% Train model
%--------------------------------------------------------------------------

SegModel('train',dir_data,opt);