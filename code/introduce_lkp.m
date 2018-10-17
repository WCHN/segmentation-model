function varargout = introduce_lkp(dat,model,opt,varargin)

K  = opt.template.K;
S0 = numel(dat);   

if S0 == 1
    GaussPrior = varargin{1};
    
    modality = dat.modality{1}.name; 
    lkp      = get_par('lkp',modality,opt);
        
    if numel(lkp) == K
        varargout{1} = dat;
        varargout{2} = model.GaussPrior(dat.population);
    
        return
    end
    
    % Modify posteriors
    gmm{1} = dat.gmm.cluster{1}{1};
    gmm{2} = dat.gmm.cluster{1}{2};
    gmm{3} = dat.gmm.cluster{2}{1};
    gmm{4} = dat.gmm.cluster{2}{2};

    [gmm,mg] = spm_gmm_lib('extras', 'more_gmms', gmm, lkp);           

    dat.gmm.cluster{1}{1} = gmm{1};
    dat.gmm.cluster{1}{2} = gmm{2};
    dat.gmm.cluster{2}{1} = gmm{3};
    dat.gmm.cluster{2}{2} = gmm{4};

    dat.gmm.part.lkp = lkp;
    dat.gmm.part.mg  = mg;
    
    % Modify GaussPrior   
    GaussPrior{7}         = lkp;
    GaussPrior{6}.ElnDetV = zeros(1,numel(lkp));                
    GaussPrior(1:4)       = spm_gmm_lib('extras', 'more_gmms', GaussPrior(1:4), lkp);        
    
    varargout{1} = dat;
    varargout{2} = GaussPrior;
else
    it_mod = varargin{1};
    
    if it_mod == opt.start_it.do_mg
        populations  = spm_json_manager('get_populations',dat);
        P            = numel(populations);

        for p=1:P  
            population0 = populations{p}.name;
            modality    = populations{p}.type;
            GaussPrior  = model.GaussPrior(population0);
            lkp         = get_par('lkp',modality,opt);

            if numel(lkp) == K
                continue
            end

            mg               = ones(1,numel(lkp));
            for k=1:max(lkp)
                kk           = sum(lkp == k);
                mg(lkp == k) = 1/kk;
            end   

            % Modify subject-specific posteriors
            for s=1:S0
                population = dat{s}.population;

                if strcmp(population0,population)
                    dat{s}.gmm.part.lkp = lkp;
                    dat{s}.gmm.part.mg  = mg;

                    gmm{1} = dat{s}.gmm.cluster{1}{1};
                    gmm{2} = dat{s}.gmm.cluster{1}{2};
                    gmm{3} = dat{s}.gmm.cluster{2}{1};
                    gmm{4} = dat{s}.gmm.cluster{2}{2};

                    gmm = spm_gmm_lib('extras', 'more_gmms', gmm, lkp);           

                    dat{s}.gmm.cluster{1}{1} = gmm{1};
                    dat{s}.gmm.cluster{1}{2} = gmm{2};
                    dat{s}.gmm.cluster{2}{1} = gmm{3};
                    dat{s}.gmm.cluster{2}{2} = gmm{4};
                end
            end

            % Modify GaussPrior   
            GaussPrior{7}         = lkp;
            GaussPrior{6}.ElnDetV = zeros(1,numel(lkp));                                        
            GaussPrior(1:4)       = spm_gmm_lib('extras', 'more_gmms', GaussPrior(1:4), lkp);           

            model.GaussPrior(population0) = GaussPrior;
        end
    end
    
    varargout{1} = dat;
    varargout{2} = model;
end
%==========================================================================