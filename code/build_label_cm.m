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
    
    %----------------------------------------------------------------------
    % Build confusion matrix 
    %----------------------------------------------------------------------
    
    ix_ul = max(ix) + 1;    % Unlabelled voxels (0) are assumed to have value: num(labels) + 1
    CM    = zeros(ix_ul,K); % Confusion matrix is of size: (num(labels) + 1) x num(tissue)
    
    % For the unlabelled voxels, we assume uniform probability
    ix_k            = ix == 0; 
    CM(ix_ul, ix_k) = rater_sens/nnz(ix_k); 
    CM(ix_ul,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
    
%     CM(ix_ul, :) = 1; 
%     prob_bg         = 0.8;
%     ix_k            = ix == 0;
%     CM(ix_bg, ix_k) = rater_sens/nnz(ix_k);
%     CM(ix_bg,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
%     CM(ix_bg, ix_k) = prob_bg;
%     CM(ix_bg,~ix_k) = 1 - prob_bg;    

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
%     CM = bsxfun(@rdivide,CM,sum(CM,2));
    
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