function [dat,model,opt] = init_all(dat,opt)

[dat,opt]   = build_dir_structure(dat,opt); % Build directory structure
[dat,opt]   = init_reg(dat,opt);            % Init registration
dat         = init_bf(dat,opt);             % Init bias field parameters
[dat,~,opt] = build_label_cm(dat,opt);      % Build label confusion matrix
[~,dat]     = init_objval(dat);             % Init objective value structs, for tracking lower bounds
dat         = init_armijo(dat);             % Init Gauss-Newton step-size, which will change depending on convergence
dat         = init_lkp(dat,opt);            % Init as one Gaussian per tissue
dat         = init_mrf(dat,opt);            % Init use of first-order MRF
% Init template and GMM
if opt.template.do
    % Init template from histogram representations of input images
    model       = init_uniform_template(dat,opt); % Create initial uniform template     
    [dat,model] = init_gmm(dat,model,opt);
    [model,dat] = update_template(dat,model,opt,true);
    
    if opt.verbose.model >= 3, show_segmentations(dat,opt); end
    if opt.verbose.model >= 3, show_PropPrior(dat,model,opt); end
else  
    % When segmenting a single subject
    [dat,model,opt] = load_model(dat,opt); % Get model parameters (model)        
end
%==========================================================================