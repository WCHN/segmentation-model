function show_template(model,opt,S0,figname)

if nargin < 4
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

title_name = 'Template';
if nargin >= 3
    title_name = [title_name ' (S=' num2str(S0) ')'];
end

spm_gmm_lib('plot','showcatimg',mu,{title_name},opt.model.nam_cls);

deal_figs(model);
%==========================================================================