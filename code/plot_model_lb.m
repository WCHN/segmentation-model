function plot_model_lb(dat,model,it,opt)

% Get figure (create if it does not exist)
figname = '(SPM) Model lower bound';
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

P = 0;
if opt.gmm.GaussPrior.constrained
    populations = spm_json_manager('get_populations',dat);
    P           = numel(populations);
end

% Plot
nfigs = 4 + P;
nrows = floor(sqrt(nfigs));
ncols = ceil(nfigs/nrows);             

subplot(nrows,ncols,1)
plot(model.lb(1:end))
title(['Model (iter=' num2str(it) ')'])

subplot(nrows,ncols,2)
plot(model.template.objval.post(1:end))
title('-ln(p(a|.))')

subplot(nrows,ncols,3)
plot(model.template.objval.likel(1:end))
title('-ln(p(.|a))')

subplot(nrows,ncols,4)
plot(model.template.objval.pr(1:end))
title('-ln(p(a))')

if opt.gmm.GaussPrior.constrained
    for p=1:P
        name = populations{p}.name;
        pr   = model.GaussPrior(name);

        subplot(nrows,ncols,4 + p)
        plot(pr{6}.KL_qVpV(1:end))
        title(['KL(qV|pV) (' name ')'])
    end
end

drawnow;

deal_figs(model);
%==========================================================================