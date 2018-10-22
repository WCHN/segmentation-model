% This script demonstrates using the segmentation-model code to learn from
% MR scans part of the MRBrainS18 segmentation challenge:
%     http://mrbrains18.isi.uu.nl/
%
% The code learns a template, hyper-parameters of a Gaussian mxiture model
% for intensities, and hyper-parameters of a dirichlet prior on the
% proportions of tissues. All these are stored in './output/train/model'.
%
% Once trained, the model can be used to segment a new subject using the
% seg_MRBrainS18_2d script in this folder.
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

% Directory for learnt model (and temp files)
opt.dir_output = './output';

% Template options
opt.template.do       = true;
opt.model.nam_cls     = {'1.gm','2.bas','3.wm','4.wmh','5.csf','6.ven','7.bg','8.st','9.skl'};
opt.template.K        = 9;
opt.template.bg_class = 7;

% For mapping labels to tissue classes
map                = containers.Map;
map('MRBrainS18')  = [1 2 3 4 5 6 0 0 0];
opt.gmm.labels.cm  = map;
opt.gmm.labels.use = true;

% For using multiple Gaussians per tissue
map          = containers.Map;
map('MRI')   = [1 2 3 4 5 5 6 6 7 7 8 8 8 9 9];
opt.dict.lkp = map;

%--------------------------------------------------------------------------
% Train model
%--------------------------------------------------------------------------

SegModel('train',dir_data,opt);