function model = update_PropPrior(dat,model,opt,it_mod)
% FORMAT [alpha,gn] = update_dirichlet_prior(alpha, meanLogX)
% 
% alpha    - Kx1 - Previous value of the Dirichlet parameters
% meanLogX - Kx1 - Mean of the log of the observations (mean(log(X)))

if ~opt.model.PropPrior.do
    return;
end

if it_mod < opt.start_it.do_prop
    return;
end

S0    = numel(dat);
alpha = model.PropPrior.alpha;
K     = numel(alpha);

meanLogX = 0;
for s=1:S0
    meanLogX = meanLogX + log(spm_matcomp('softmax',dat{s}.gmm.prop));
end
meanLogX = meanLogX./S0;

alpha    = double(alpha(:));
logalpha = log(alpha);
meanLogX = double(meanLogX(:));
K        = numel(alpha);

E  = NaN;
for gn=1:100000
    
    % Compute objective function (negtive log-likelihood)
    Eo = E;
    E  = sum(gammaln(alpha)) ...
         - gammaln(sum(alpha)) ...
         - sum((alpha-1).*meanLogX);
%     fprintf('E = %10.6g\n', E);
    if E < Eo && abs(Eo - E) < 1E-7
        % It sometimes overshoots during the first iterations but gets back
        % on track on its own.
        break
    end
    
    % Compute grad/hess
    g = alpha .* ( psi(alpha) - psi(sum(alpha)) - meanLogX );
    H = (alpha * alpha') .* (diag(psi(1,alpha)) - psi(1,sum(alpha)) * ones(K));
    
    H = spm_matcomp('LoadDiag', H);
    
    % update
    logalpha = logalpha - H\g;
    alpha    = exp(logalpha);
end

model.PropPrior.alpha = alpha;
model.PropPrior.norm  = S0*(gammaln(sum(alpha)) - sum(gammaln(alpha)));

% Show results
show_PropPrior(dat,model,opt);

% Save updated PropPrior
PropPrior = model.PropPrior;
fname     = fullfile(opt.dir_model,'PropPrior.mat');
save(fname,'PropPrior');
%==========================================================================