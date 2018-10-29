function [dat,model,scl] = init_gmm(dat,model,opt)

% Parameters
%--------------------------------------------------------------------------
S0         = numel(dat);
niter      = opt.gmm.hist.niter_main;
K          = opt.template.K;
verbose    = opt.gmm.hist.verbose;
nii_a      = model.template.nii;
tiny       = get_par('tiny');
opt.gmm.verbose = opt.gmm.hist.verbose_gmm;

% Observations from histograms
%--------------------------------------------------------------------------
[dat,obs,scl,bw] = get_hists(dat,opt);

% Init VB-GMM
%--------------------------------------------------------------------------
[dat,model] = get_gmms(obs,model,dat,opt);

% Learn Gauss-Wishart hyper-parameters
%--------------------------------------------------------------------------
lb = -Inf(S0,niter);
for iter=1:niter
    		
%     for s=1:S0
    parfor s=1:S0 % Iterate over subjects
        population = dat{s}.population;          
        
        % Fit VBGMM (update posteriors and mixing weights)
        %------------------------------------------------------------------  
        pr                  = model.GaussPrior(population);
        [dat{s},lb(s,iter)] = update_gmm_hist(obs{s},dat{s},bw{s},pr,opt);                                       
    end    
    
	% Update Gauss-Wishart hyper-parameters
	%----------------------------------------------------------------------    
    model = update_GaussPrior(dat,model,opt);    
end

if opt.template.load_a_der            
    % Compute template derivatives (will be loaded later in
    % update_template)
    %----------------------------------------------------------------------
        
    dat   = init_load_a_der(dat,opt);
    mat_a = nii_a.mat;                     
    dm_a  = nii_a.dat.dim;  
      
    parfor s=1:S0
%     for s=1:S0                
        [obs,dm_s,mat_s,vs_s,scl1] = get_obs(dat{s},'do_scl',true);                      
        labels                     = get_labels(dat{s},opt);                                
        miss                       = get_par('missing_struct',obs);
        % figure; imshow3D(reshape(labels{1},dm_s))

        % Set proportions
        prop            = dat{s}.gmm.prop;         
        dat{s}.gmm.prop = log(prop); % Make logarithmic
        
        % Get responsibilities      
        Z = get_resp(obs,1,dat{s},[],labels,scl1,miss,dm_s,opt);
        labels = []; obs = [];
        
        if 0
            figure(666); k = 1; imshow3D(reshape(Z(:,:,:,k),[dm_s(1:2) dm_s(3)]))   
            figure(666); imshow3D(squeeze(reshape(Z,[dm_s K]))) 
        end
        
        % Compute warp from template to subject space
        E      = spm_dexpm(dat{s}.reg.r,opt.reg.B);
        Affine = mat_a\E*mat_s;                
        y      = spm_warps('transform',Affine,spm_warps('identity',dm_s));
        if dm_s(3) == 1
            y(:,:,:,3) = 1;
        end

        % Push responsibilities in subject space to template space
        [Z,dat{s}.template.bb] = push_responsibilities(Z,y,dm_a(1:3),mat_a,dat{s}.template.bb);
        y                      = [];        
            
        % Compute gradients and Hessian
        a          = rotate_template(nii_a,opt);
        [gr,H,dll] = diff_template(a,Z,dat{s}.gmm.prop,dat{s}.template.bb,vs_s,opt); 
        Z = []; a = [];   

        dat{s}.template.ll = dll;                    
        
        % Write derivatives to disk                                    
        spm_misc('create_nii',dat{s}.template.pth_gr,gr,mat_a,[spm_type('float32') spm_platform('bigend')],'a_gr'); gr = [];
        spm_misc('create_nii',dat{s}.template.pth_H,H,mat_a,[spm_type('float32') spm_platform('bigend')],'a_H');    H  = [];
    end
end

% Set ElnDetV to zero
populations  = spm_json_manager('get_populations',dat);
P            = numel(populations);
for p=1:P  
    population = populations{p}.name;
    modality   = populations{p}.type;
    GaussPrior = model.GaussPrior(population);
    lkp        = get_par('lkp',modality,opt);
    
    GaussPrior{6}.ElnDetV = zeros(1,max(lkp));                
    
    model.GaussPrior(population) = GaussPrior;
end
%==========================================================================
    
%==========================================================================                            
function [dat,lb] = update_gmm_hist(obs,dat,scl,GaussPrior,opt)
Verbose = opt.gmm.verbose;
cluster = dat.gmm.cluster;
prop    = dat.gmm.prop; 
IterMax = opt.gmm.hist.niter_gmm;
C       = size(obs{1},2);
prop    = {'Prop',prop};
 
BinUncertainty = ((scl.*ones(1,C)).^2)./12; % Variance of uniform distribution

[~,clust,prp,lb] = spm_gmm_loop(obs,cluster,prop, ...
                                'GaussPrior',GaussPrior, ...
                                'BinUncertainty',BinUncertainty, ...
                                'IterMax',IterMax,'SubIterMax',IterMax, ...
                                'Verbose',Verbose);

lb              = lb.sum(end);
dat.gmm.cluster = {{clust.MU,clust.b},{clust.V,clust.n}};
dat.gmm.prop    = prp.Prop; 
%==========================================================================

%==========================================================================
function [dat,obs,scl,bw] = get_hists(dat,opt)
S0  = numel(dat);
obs = cell(S0,1);

% Fit histograms
%--------------------------------------------------------------------------
scl = cell(1,S0);
bw  = cell(1,S0);
for s=1:S0  
    modality             = dat{s}.modality{1}.name;
    [obs_s,~,~,~,scl{s}] = get_obs(dat{s},'do_scl',true);
    obs_s                = double(obs_s);
    C                    = size(obs_s,2);
    
    % Fit histogram to image intensities
    %----------------------------------------------------------------------
    mn = min(obs_s,[],1);
    mx = max(obs_s,[],1);
    if strcmpi(modality,'ct') 
        bins  = single((mn:mx).');
        bw{s} = mean(diff(bins));        
    else
        bins  = cell(1,C);
        nbins = 64;%floor(1000/2^(1 + C));
        for c=1:C
            bins{c}  = single(linspace(mn(c),mx(c),nbins)');
            bins{c}  = round(bins{c});
            bw{s}(c) = mean(diff(bins{c}));
        end        
    end
            
    labels = get_labels(dat{s},opt);
    
    [V,W] = spm_imbasics('hist',obs_s,bins,'KeepZero',false,'Missing',true,'Labels',labels);
        
    obs{s} = {single(V),single(W)};
end
%==========================================================================    