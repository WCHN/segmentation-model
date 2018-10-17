function [dat,CM,opt] = build_label_cm(dat,opt)

% Parameters
cm         = opt.gmm.labels.cm;
rater_sens = opt.gmm.labels.S;
S0         = numel(dat);
K          = opt.template.K;
CM         = 0;

% For each subject with labels, build the label confusion matrix
for s=1:S0
    population = dat{s}.population;  
    
    if ~cm.isKey(population)
        dat{s}.gmm.cm = zeros(1,K);
        continue;
    else
        ix = cm(population);
    end
    
    if numel(ix)~=K
        error('numel(cm_p)~=K')
    end
    
    % Build confusion matrix (max(ix) x K)    
    ix_bg = max(ix) + 1;
    CM    = zeros(ix_bg,K);    
    
    % Background labels (assumed to always have the label 0)
    ix_k            = ix == 0;
    CM(ix_bg, ix_k) = rater_sens/nnz(ix_k);
    CM(ix_bg,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
    
    % Other labels
    for k=ix
        if k == 0, continue; end
        
        ix_k        = ix == k;
        CM(k, ix_k) = rater_sens/nnz(ix_k);
        CM(k,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
    end

    CM = bsxfun(@rdivide,CM,sum(CM,1));
    
%     % Adjust CM for when using multiple Gaussians per tissue
%     modality = dat{s}.modality{1}.name;
%     if opt.dict.lkp.isKey(modality)
%         lkp = opt.dict.lkp(modality);    
%     else
%         lkp = 1:K;
%     end
%     
%     w = zeros(1,K);
%     for k=1:K
%         w(k) = 1/sum(lkp == k);
%     end
%     
%     CM = bsxfun(@times,CM,w);
    
    % Assign to dat
    dat{s}.gmm.cm = CM;
end
%===========================================================================