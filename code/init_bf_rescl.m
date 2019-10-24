function dat = init_bf_rescl(dat,opt,model)
% FORMAT dat = init_bf_rescl(dat,opt,model)
% dat   - Subjects data structure
% opt   - Options structure
%
% Initialise mapping between GMM clusters and Template classes:
% * dat.gmm.part.lkp: GMM to Template mapping
% * dat.gmm.part.mg:  Within class weighting
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Create scaling for each population       
populations = spm_json_manager('get_populations',dat);    
P           = numel(populations);
S0          = numel(dat);
map         = containers.Map;
val_def     = 100;

if nargin == 3
    % If segmenting, not learning model, make sure that the bias field
    % scaling is based on the intensity prior.

    if ~iscell(model)
        model = {model};
    end
    
    % Sum over template voxels
    template = single(model{1}.template.nii.dat(:,:,:,:));
    template = spm_matcomp('softmax',template);  
    dm_a     = size(template);
    K        = dm_a(4);
    template = reshape(template,[prod(dm_a(1:3)) K]);
    sm       = sum(template,1);
    sm       = sm./sum(sm);   
    clear template

    for p=1:P % Loop over populations

        % Get population intensity prior
        nam_population = populations{p}.name;
        modality       = get_modality_name(dat,nam_population);
        C              = populations{p}.C;
        if strcmpi(modality,'CT') 
            val_p = ones(C,1); 
        else
            pr  = model{1}.GaussPrior(nam_population);
            MU0 = pr{1};
            lkp = pr{7};

            % Create mean class weights and rescale means
            mg = ones(1,numel(lkp));
            for k=1:max(lkp)
                kk           = sum(lkp == k);
                mg(lkp == k) = 1/kk;
            end   
            MU0 = mg.*MU0;

            % Map many-means-per-tissue to one-mean-per tissue
            val_p = zeros(C,K);
            for k=1:K
                val_p(:,k) = sum(MU0(:,lkp == k),2);
            end

            % Create scaling
            val_p = sm.*val_p;    
            val_p = sum(val_p,2);            
        end
        
        map(nam_population) = val_p;        
    end        
elseif nargin < 3
    
    for p=1:P % Loop over populations
                
        nam_population = populations{p}.name;
        modality       = get_modality_name(dat,nam_population);
        C              = populations{p}.C;
        if strcmpi(modality,'CT') 
            val_p = ones(C,1);        
        else
            val_p = val_def*ones(C,1);
        end
        
        map(nam_population) = val_p;        
    end
end

% Add scaling to each subjects dat struct
for p=1:P
    nam_pop = populations{p}.name;        

    for s=1:S0
        nam_pop_s = dat{s}.population;

        if strcmp(nam_pop,nam_pop_s)
            dat{s}.bf.scl = map(nam_pop_s);
        end
    end
end
%==========================================================================