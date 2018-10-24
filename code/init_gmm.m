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
[dat,model] = get_gmms(obs,model,dat,K,opt);

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
function [dat,model] = get_gmms(obs,model,dat,K,opt)
S0      = numel(obs);
init_ix = opt.gmm.hist.init_ix;
tiny    = get_par('tiny');
        
% Set posteriors
%--------------------------------------------------------------------------
for s=1:S0    
    modality    = dat{s}.modality{1}.name;   
    dm          = obs_info(dat{s}); 
    lkp         = get_par('lkp',modality,opt);
    [~,ix_tiny] = get_par('ix_tiny',dat{s}.population,lkp,opt);
    
    % Initialise mixing proportions    
    prop            = ones(1,K)./K;   
    dat{s}.gmm.prop = prop;        

    if strcmpi(modality,'ct')
        
        %------------------------------------------------------------------
        % Initilisation of GMMs when data is CT
        %------------------------------------------------------------------                

        SHOW_FIT = 0;
        
        X  = double(obs{s}{1});
        B  = double(obs{s}{2});

        [~,MU,~,~,b,V,n] = spm_gmm(X,K,B,'BinWidth',1,'GaussPrior',opt.ct.GaussPrior,'PropPrior',dat{s}.gmm.prop,'Start','prior','Verbose',SHOW_FIT,'IterMax',10);
        
        post = {{MU,b},{V,n}}; % Posteriors        
    else        
        
        %------------------------------------------------------------------
        % Initilisation of GMMs when data is MRI
        %------------------------------------------------------------------
        
        X = get_obs(dat{s},'do_scl',true);    
        C = size(X,2);
        
        miss.C  = spm_gmm_lib('obs2code',X);
        miss.L  = unique(miss.C);
        miss.nL = numel(miss.L);

        if opt.gmm.labels.use && isfield(dat{s},'label') && opt.gmm.labels.cm.isKey(dat{s}.population)                        
            
            %--------------------------------------------------------------
            % Labels provided
            %--------------------------------------------------------------        
            
            sort_pars = false;
            
            ix = opt.gmm.labels.cm(dat{s}.population);
                                    
            % Get labels
            %--------------------------------------------------------------
            labels = get_labels(dat{s},opt);                        
            labs   = ix(ix~=0);
            K_l    = numel(labs); % Number of labelled classes
            K_nl   = K - K_l;     % Number of non-labelled classes
            ix_bg  = max(ix) + 1;
                
            % Get resps for non-labelled (background) voxels
            msk_nl = labels{1} == ix_bg;
            % figure; imshow3D(squeeze(reshape(msk_nl ,[dm_s])))
            X1     = X(msk_nl,:);
            
            Z_nl = resp_from_kmeans(X1,K_nl - sum(ix_tiny),'uniform');                        
            clear X1
            
            % Get resps for labelled voxels
            cnt = 1;
            Z   = zeros([size(X,1) K],'single');
            for k=1:K
                if ix_tiny(k)
                    continue; 
                end
                
                if ix(k)==0
                    Z(msk_nl,k) = Z_nl(:,cnt);
                    cnt         = cnt + 1;
                else    
                    Z(:,k)      = labels{1}==ix(k);
                end
            end
            clear Z_nl msk_nl                                                
        else                      
            %--------------------------------------------------------------
            % No labels provided
            % Get responsibilities from kmeans labels
            %--------------------------------------------------------------
            
            sort_pars = true;
            
            Z = resp_from_kmeans(X,K - sum(ix_tiny),'plus');
                        
            nZ  = zeros([size(X,1) K],'single');
            cnt = 1;
            for k=1:K
                if ix_tiny(k)
                    continue; 
                end
                
                nZ(:,k) = Z(:,cnt);
                cnt     = cnt + 1;
            end
            Z = nZ;
            clear nZ
        end              
        
        % Add tiny value where resps are zero
        Z(Z==0) = tiny;
        Z       = bsxfun(@rdivide,Z,sum(Z,2));                    
        % figure; imshow3D(squeeze(reshape(Z(:,:),[dm K]))) 

        [SS0,SS1] = spm_gmm_lib('SuffStat',X,Z,1);
        
        % Build prior
        b  = ones(1,K);
        n  = C*ones(1,K);
        MU = double(bsxfun(@rdivide,SS1,SS0));                
        V  = C*eye(C);
        V  = repmat(V,[1 1 K]);
        
        % Set prior
        pr = {MU,b,V,n};

        % Compute suffstats from responsibilities, then GMM parameters
        post = get_cluster(X,ones(size(X),'single'),dm,pr,miss,Z,'sort_pars',sort_pars);                
    end       

    dat{s}.gmm.cluster = post; % Posteriors
end

% Set intensity prior as sample mean (for MRI)
%--------------------------------------------------------------------------
populations   = spm_json_manager('get_populations',dat);
P             = numel(populations);
GaussPrior    = containers.Map;
lb_pr         = struct;                                        
lb_pr.KL_qVpV = 0;
lb_pr.ElnDetV = zeros(1,K);
for p=1:P  
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    names       = get_channel_names(dat,populations{p});
        
    if strcmpi(modality,'CT')
        GaussPrior(population0) = opt.ct.GaussPrior;
    else
        % Use sample mean
        MU  = 0;
        b   = 0;
        W   = 0;
        n   = 0;
        cnt = 0;
        for s=1:S0
            population = dat{s}.population;

            if strcmp(population0,population)
                MU = MU + dat{s}.gmm.cluster{1}{1};
                b  = b  + dat{s}.gmm.cluster{1}{2};
                W  = W  + dat{s}.gmm.cluster{2}{1};
                n  = n  + dat{s}.gmm.cluster{2}{2};

                cnt = cnt + 1;
            end
        end

        MU = MU./cnt;
        b  = b./cnt;
        W  = W./cnt;
        n  = n./cnt;
        
        lkp = 1:K;
        
        GaussPrior(population0) = {MU,b,W,n,names,lb_pr,lkp};            
    end
end
model.GaussPrior = GaussPrior;
clear GaussPrior

% Compute initial estimate of hyper-parameters of VB-GMM
%--------------------------------------------------------------------------
model = update_GaussPrior(dat,model,opt);

% Set tissue indices to match between different populations
%--------------------------------------------------------------------------
for p=1:P
    population0 = populations{p}.name;  
        
    if init_ix.isKey(population0) 
        ix = init_ix(population0);    
    else                          
        continue;
    end
    
    % Adjust posteriors
    for s=1:S0
        population = dat{s}.population;
        
        if strcmp(population0,population)
            dat{s}.gmm.cluster{1}{1} = dat{s}.gmm.cluster{1}{1}(:,ix);
            dat{s}.gmm.cluster{1}{2} = dat{s}.gmm.cluster{1}{2}(ix);
            dat{s}.gmm.cluster{2}{1} = dat{s}.gmm.cluster{2}{1}(:,:,ix);
            dat{s}.gmm.cluster{2}{2} = dat{s}.gmm.cluster{2}{2}(ix);            
        end
    end
    
    % Adjust priors
    pr    = model.GaussPrior(population0);
    pr{1} = pr{1}(:,ix);
    pr{2} = pr{2}(ix);
    pr{3} = pr{3}(:,:,ix);
    pr{4} = pr{4}(ix);
    model.GaussPrior(population0) = pr;    
end
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
    type                 = dat{s}.modality{1}.name;
    [obs_s,~,~,~,scl{s}] = get_obs(dat{s},'do_scl',true);
    obs_s                = double(obs_s);
    C                    = size(obs_s,2);
    
    % Fit histogram to image intensities
    %----------------------------------------------------------------------
    mn = min(obs_s,[],1);
    mx = max(obs_s,[],1);
    if strcmpi(type,'ct') 
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
        
    obs{s} = {V,W};
end
%==========================================================================
     
%==========================================================================   
function Z = resp_from_kmeans(X,K,start_method)
% Get initial labels using kmeans
L = spm_kmeans(X,K,'Distance','cityblock','Start',start_method);

% Compute responsibilities from kmeans labels
Z = zeros([numel(L) K],'single');
for k=1:K
    Z(:,k) = L(:)==k; 
end
%==========================================================================   

%==========================================================================   
function chn_names = get_channel_names(dat,population)
S0    = numel(dat);
C     = population.C; 
name  = population.name; 
for s=1:S0
   population0 = dat{s}.population;
   
   if strcmpi(population0,name)              
       [~,~,~,~,~,~,~,chn_names] = obs_info(dat{s});              
       
       return
   end
end
%==========================================================================       
