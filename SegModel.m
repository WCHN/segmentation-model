function varargout = SegModel(varargin)
%__________________________________________________________________________
%
% Code for segmentation model
%
%--------------------------------------------------------------------------
%
% SegModel('train',dat,opt)
% > Train segmentation model from a bunch of images
%
% res = SegModel('segment',dat,opt)
% > Segment images with trained segmentation model
%
% Arguments
% ----------
% dat - [1xS] Cell array with subject data
% opt -       Model options
% res - [1xS] Cell array with segmentation results
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin == 0
    help run_SegModel
    error('Not enough argument. Type ''help spm_gmm_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'train'
        [varargout{1:nargout}] = SegModel_train(varargin{:});        
    case 'segment'
        [varargout{1:nargout}] = SegModel_segment(varargin{:});           
    otherwise
        help run_SegModel
        error('Unknown function %s. Type ''help run_SegModel'' for help.', id)
end
%==========================================================================

%==========================================================================
function SegModel_train(dat,opt)
% Train segmentation model from a bunch of images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

dat = SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Parse options
%--------------------------------------------------------------------------

[opt,holly] = default_opt(opt);

%--------------------------------------------------------------------------
% Initialise model
%--------------------------------------------------------------------------

[dat,model,opt] = init_all(dat,opt);

%--------------------------------------------------------------------------
% Train model
%--------------------------------------------------------------------------

model.lb = -Inf;
for it_mod=1:opt.model.niter                            
    
    % Set-up loading template derivatives from disk       
    dat = init_load_a_der(dat,opt);

    % Introduce using multiple Gaussians per tissue
    [dat,model] = introduce_lkp(dat,model,opt,it_mod);

    % Do segmentation
    [holly,dat] = distribute(holly,'segment_subject','inplace',dat,model,opt,it_mod);   

    % Save dat (for debugging)
    fname = fullfile(opt.dir_model,'dat.mat');
    save(fname,'dat','-v7.3');

    % Some model parameters changes with iteration number
    opt = modify_opt(opt,it_mod);                      

    % Mean correct the rigid-body transforms
    dat = meancorrect_aff(dat,opt);

    % Mean correct bias field (also updates posteriors)
    dat = meancorrect_bf(dat,model.GaussPrior,opt,it_mod);          

    % Update template	    
    model = update_template(dat,model,opt);           

    % Update Gauss-Wishart hyper-parameters	    
    model = update_GaussPrior(dat,model,opt,it_mod); 

    % Update proportions hyper-parameter
    model = update_PropPrior(dat,model,opt,it_mod);

    % Compute model lower bound
    model = compute_model_lb(dat,model,opt);     

    % Check convergence
    gain = spm_misc('get_gain',model.lb);
    if opt.verbose.model >= 1
        fprintf('%3s | %3d | lb = %10.6f | gain = %1.5f\n', 'mod', it_mod, model.lb(end), gain);
    end        

    % Some verbose
    if opt.verbose.model >= 2, plot_model_lb(dat,model,it_mod,opt); end
    if opt.verbose.model >= 3, show_segmentations(dat); end
end

if opt.model.clean_up
    % Clean-up temporary files
    rmdir(opt.dir_vel,'s');
    rmdir(opt.dir_a_der,'s');
    rmdir(opt.dir_seg2d,'s');
end
%==========================================================================

%==========================================================================
function res = SegModel_segment(dat,opt)
% Segment images with trained segmentation model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

dat = SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Parse options
%--------------------------------------------------------------------------

[opt,holly] = default_opt(opt);

%--------------------------------------------------------------------------
% Initialise model
%--------------------------------------------------------------------------

[dat,model,opt] = init_all(dat,opt);

%--------------------------------------------------------------------------
% Segment subject(s)
%--------------------------------------------------------------------------
               
[~,~,res] = distribute(holly,'segment_subject','inplace',dat,'iter',model,opt);  
%==========================================================================

%==========================================================================
function dat = SegModel_init(dat,opt)
% Initialise segmentation model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if ~isfield(opt.dep,'aux_toolbox')
    error('~isfield(opt.dep,''aux_toolbox'')')
else
    if ~isdeployed, addpath(opt.dep.aux_toolbox);  end
end

if ~isfield(opt.dep,'dist_toolbox')
    error('~isfield(opt.dep,''dist_toolbox'')')
else
    if ~isdeployed, addpath(opt.dep.dist_toolbox); end
end

if isfield(opt,'clean') && isfield(opt.clean,'cnn_mrf') && isfield(opt.clean.cnn_mrf,'do')
    if opt.clean.cnn_mrf.do
        if ~isfield(opt.dep,'cnn_mrf')
            error('~isfield(opt.dep,''cnn_mrf'')')
        else
            if ~isdeployed, addpath(opt.dep.cnn_mrf); end
        end
    end
end

if ~isdeployed, addpath(fullfile(fileparts(mfilename('fullpath')),'code')); end

set_globals;

if ischar(dat)
    dat = spm_json_manager('init_dat',dat);
    if ~(isfield(opt,'template') && isfield(opt.template,'do') && opt.template.do)
        dat = dat(1);
    end
end
%==========================================================================