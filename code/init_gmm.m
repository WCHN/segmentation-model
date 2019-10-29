function [dat,model,scl] = init_gmm(dat,model,opt)
% FORMAT [dat,model,scl] = init_gmm(dat,model,opt)
% dat   - Subjects data structure
% model - Model structure
% opt   - Options structure
% scl   - 
%
% Initialise GMM parameters by fitting a GMM to the image histograms.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
%--------------------------------------------------------------------------
S0                             = numel(dat);
niter                          = opt.gmm.hist.niter_main;
opt.gmm.verbose                = opt.gmm.hist.verbose_gmm;
opt.gmm.GaussPrior.constrained = false; % This is disabled here just for speed

if ~opt.given.GaussPrior && ~opt.GaussPrior.uninformative
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
else    
    K            = opt.template.K;
    lkp0         = 1:K;
    populations  = spm_json_manager('get_populations',dat);
    P            = numel(populations);
        
    if opt.GaussPrior.uninformative
        
        %----------------
        % Init uninformative GaussPrior
        %----------------
        
        GaussPrior = containers.Map;
        
        for p=1:P  
            population = populations{p}.name;
            C          = populations{p}.C;
            names      = get_channel_names(dat,populations{p});
            
            b  = ones(1,K);
            n  = C*ones(1,K);
            MU = zeros(C,K);        
            V  = repmat(eye(C),[1 1 K]);

            pr{1} = MU;
            pr{2} = b;
            pr{3} = V;
            pr{4} = n;    
            pr{5} = names;    
            lb_pr         = struct;                                        
            lb_pr.KL_qVpV = 0;
            lb_pr.ElnDetV = zeros(1,K);
            pr{6}         = lb_pr;
            pr{7} = lkp0;
        
            GaussPrior(population) = pr;
        end
                                
    else
        
        %----------------
        % Load GaussPrior
        %----------------
        
        load(opt.gmm.pth_GaussPrior)

        % Collapse prior        
        for p=1:P  
            population = populations{p}.name;
            pr         = GaussPrior(population);
            lkp        = pr{7};

            if ~isequal(lkp0,lkp) 
                pr(1:4)       = spm_gmm_lib('extras','collapse_gmms',pr(1:4),lkp);
                lb_pr         = struct;                                        
                lb_pr.KL_qVpV = 0;
                lb_pr.ElnDetV = zeros(1,numel(lkp0));
                pr{6}         = lb_pr;
                pr{7}         = lkp0;
            end        

            GaussPrior(population) = pr;
        end       
    end
    
    model.GaussPrior = GaussPrior;
end

% Set prop to uniform
for s=1:S0
    dat{s}.gmm.prop = ones(1,opt.template.K)./opt.template.K;
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
parfor s=1:S0  
% for s=1:S0, fprintf('obs! for s=1:S0\n')
    modality             = dat{s}.modality{1}.name;
    [obs_s,~,~,~,scl{s}] = get_obs(dat{s},'do_scl',true,'mskonlynan',opt.seg.mskonlynan,'missmod',opt.seg.missmod,'no_scl',opt.bf.no_scl); % Do not subsample!
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
            
%     labels = get_labels(dat{s},opt);
    
    [V,W] = spm_imbasics('hist',obs_s,bins,'KeepZero',false,'Missing',true);
        
    obs{s} = {single(V),single(W)};
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