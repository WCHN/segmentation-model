function show_template(model,opt,figname)

if nargin < 3
    figname = '(SPM) Template';
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

mu = single(model.template.nii.dat(:,:,:,:));
mu = spm_matcomp('softmax',mu);

spm_gmm_lib('plot','showcatimg',mu,{'Template'},opt.model.nam_cls);

deal_figs(model);
%==========================================================================