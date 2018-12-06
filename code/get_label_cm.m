function CM = get_label_cm(dat,opt)

population = dat.population;
iter       = opt.model.it;
cm         = opt.gmm.labels.cm;
% rater_sens = opt.gmm.labels.S;
rater_sens = opt.sched.labels(min(numel(opt.sched.labels),iter));
K          = opt.template.K;

if ~cm.isKey(population)
    ix = zeros(1,K);
else
    ix = cm(population);
end

if numel(ix)~=K
    error('numel(cm_p)~=K')
end

%----------------------------------------------------------------------
% Build confusion matrix 
%----------------------------------------------------------------------

ix_ul = max(ix) + 1;    % Unlabelled voxels (0) are assumed to have value: num(labels) + 1
CM    = zeros(ix_ul,K); % Confusion matrix is of size: (num(labels) + 1) x num(tissue)
ix_k  = ix == 0; 

% For the unlabelled voxels, we assume uniform probability
CM(ix_ul, ix_k) = rater_sens/nnz(ix_k); 
CM(ix_ul,~ix_k) = (1 - rater_sens)/nnz(~ix_k);

% prob_bg         = 0.5;
% CM(ix_ul, ix_k) = prob_bg;
% CM(ix_ul,~ix_k) = 1 - prob_bg;    

%     ix_k            = ix == 0;
%     CM(ix_ul, ix_k) = rater_sens/nnz(ix_k);
%     CM(ix_ul,~ix_k) = (1 - rater_sens)/nnz(~ix_k);

% For the labelled voxels, we define probabilities from a user-given
% rater sensitivity
for k=ix
    if k == 0
        % Skip unlabelled
        continue; 
    end

    ix_k        = ix == k;
    CM(k, ix_k) = rater_sens/nnz(ix_k);
    CM(k,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
end

% Normalise confusion matrix
CM = bsxfun(@rdivide,CM,sum(CM,1));
% CM = bsxfun(@rdivide,CM,sum(CM,2));

% % Adjust CM for when using multiple Gaussians per tissue
% lkp = dat.gmm.part.lkp;
% w   = zeros(1,K);
% for k=1:K
%     w(k) = 1/sum(lkp == k);
% end
% CM = bsxfun(@times,CM,w);

return