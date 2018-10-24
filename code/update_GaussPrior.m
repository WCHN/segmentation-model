function model = update_GaussPrior(dat,model,opt,do_update)
if nargin < 4, do_update = true; end

dir_model   = opt.dir_model;
verbose     = opt.gmm.hist.verbose;
constrained = opt.gmm.GaussPrior.constrained;
GaussPrior  = model.GaussPrior;

S0          = numel(dat);
populations = spm_json_manager('get_populations',dat);
P           = numel(populations);

%--------------------------------------------------------------------------
% First we update for CT data, where all populations share the same GaussPrior
%--------------------------------------------------------------------------

clusters = {};  
cnt      = 1;
    
for p=1:P % Iterate over populations
    population0 = populations{p}.name;
    modality    = populations{p}.type;
        
    if ~strcmpi(modality,'ct')
        % Skip non-CT data
        continue
    end
    
    pr = GaussPrior(population0);
    
    % Get lkp
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
            lkp = dat{s}.gmm.part.lkp;
            break
        end
    end

    for s=1:S0 % Iterate over subjects      
        population = dat{s}.population; 

        if strcmpi(population0,population)
            clusters{cnt} = dat{s}.gmm.cluster;
            cnt           = cnt + 1;
       end
    end
end

[pr(1:4),extras] = spm_gmm_lib('UpdateHyperPars',clusters,pr(1:4), ...
                               'constrained',constrained,'figname','CT','verbose',verbose,'lkp',lkp);    


lb_pr                  = pr{6};                                        
lb_pr.KL_qVpV(end + 1) = extras.lb;  % KL(q(V0) || p(V0))
lb_pr.ElnDetV          = extras.ldW; % E[ln|V0|]    
pr{6}                  = lb_pr;    

GaussPrior(population0) = pr;
    
%--------------------------------------------------------------------------
% Then we update for MR data, where all populations have their individual
% GaussPriors
%--------------------------------------------------------------------------

for p=1:P % Iterate over populations
    population0 = populations{p}.name;
    modality    = populations{p}.type;
        
    if strcmpi(modality,'ct')
        % Skip CT data
        continue
    end
    
    pr           = GaussPrior(population0);
    clusters     = {};  
    cnt          = 1;

    % Get lkp
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
            lkp = dat{s}.gmm.part.lkp;
            break
        end
    end

    for s=1:S0 % Iterate over subjects      
        population = dat{s}.population; 

        if strcmpi(population0,population)
            clusters{cnt} = dat{s}.gmm.cluster;
            cnt           = cnt + 1;
       end
    end

    [pr(1:4),extras] = spm_gmm_lib('UpdateHyperPars',clusters,pr(1:4), ...
                                   'constrained',constrained,'figname',population0,'verbose',verbose,'lkp',lkp);    


    lb_pr                  = pr{6};                                        
    lb_pr.KL_qVpV(end + 1) = extras.lb;  % KL(q(V0) || p(V0))
    lb_pr.ElnDetV          = extras.ldW; % E[ln|V0|]    
    pr{6}                  = lb_pr;    

    GaussPrior(population0) = pr;
end

fname = fullfile(dir_model,'GaussPrior.mat');
save(fname,'GaussPrior');

model.GaussPrior = GaussPrior;

deal_figs(model);
%==========================================================================