function model = compute_model_lb(dat,model,opt)
S0  = numel(dat);
dlb = 0;
for s=1:S0
    dlb = dlb + dat{s}.lb.last;
end

dlb = dlb + model.template.objval.pr(end);        

if opt.gmm.GaussPrior.constrained
    % Make sure model ll is correct when having a hierarchical prior on the
    % V hyper-parameter of the VB GMM
    populations = spm_json_manager('get_populations',dat);
    P           = numel(populations);
    for p=1:P
        name = populations{p}.name;
        pr   = model.GaussPrior(name);
        dlb  = dlb + pr{6}.KL_qVpV(end);
    end
end

dlb = dlb + model.PropPrior.norm;

model.lb(end + 1) = dlb;   

if any(~isfinite(model.lb(2:end)))
    warning('~isfinite(model.lb)');
end
%==========================================================================