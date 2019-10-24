function [dat,model,opt] = init_all(dat,opt)
% FORMAT [dat,model,opt] = init_all(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
% model - Model structure
%
% * Build directory structure
% * Init registration
% * Init bias field parameters
% * Init objective value structs, for tracking lower bounds
% * Init Gauss-Newton step-size, which will change depending on convergence
% * Init one Gaussian per tissue
% * Init use of first-order MRF
% * Init template from histogram representations of input images
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

[dat,opt] = build_dir_structure(dat,opt); % Build directory structure
[~,dat]   = init_objval(dat);             % Init objective value structs, for tracking lower bounds
[dat,opt] = init_reg(dat,opt);            % Init registration
dat       = init_armijo(dat);             % Init Gauss-Newton step-size, which will change depending on convergence
dat       = init_lkp(dat,opt);            % Init as one Gaussian per tissue
dat       = init_mrf(dat,opt);            % Init use of first-order MRF
dat       = init_isct(dat);               % Flag channel as either CT or not
% Init template and GMM
if opt.template.do
    % Init template from histogram representations of input images

    opt.given.template   = false;
    opt.given.GaussPrior = false;
    if isfield(opt.template,'pth_template') && exist(opt.template.pth_template,'file')
        opt.given.template   = true;
    end
    if isfield(opt.gmm,'pth_GaussPrior') && exist(opt.gmm.pth_GaussPrior,'file')
        opt.given.GaussPrior = true;
    end
    
    model = init_template(dat,opt); % Create initial uniform template
    
    fprintf('Start making initial estimates of GMM parameters...');
    [dat,model] = init_gmm(dat,model,opt);
    fprintf('done!\n')
    
    if opt.given.template && opt.given.GaussPrior
        dat = init_bf_rescl(dat,opt,model); % Init bf rescaling values
    else
        dat = init_bf_rescl(dat,opt);       % Init bf rescaling values
    end
    dat = init_bf(dat,opt);       % Init bias field parameters                    
    
    if ~opt.given.template
        fprintf('Start making initial estimate of template...');
        [model,dat] = update_template(dat,model,opt,true);
        fprintf('done!\n')
    end
    
    if opt.verbose.model >= 3 && opt.template.update && ~opt.given.template
        show_tissues(dat,model,opt); 
    end
%     if opt.verbose.model >= 3, show_PropPrior(dat,model,opt); end
else         
    % When segmenting a single subject
    [dat,model,opt] = load_model(dat,opt); % Get model parameters (model)        
    
    dat = init_bf(dat,opt); % Init bias field parameters
end
%==========================================================================