function [dat,res] = segment_subject(dat,model,opt,it_mod)

if nargin < 4, it_mod = 1; end

set_globals;

%--------------------------------------------------------------------------
% Initialise model parameters
%--------------------------------------------------------------------------

% Get labels
labels = get_labels(dat,opt); 
K      = opt.template.K;

% Population parameters
mat_a      = model.template.nii.mat;       % Template orientation matrix                  
dm_a       = model.template.nii.dat.dim;   % Template orientation matrix   
PropPrior  = model.PropPrior;              % Proportion prior parameter
GaussPrior = model.GaussPrior(dat.population); 

% Subject parameters
[obs,dm_s,mat_s,vs_s,scl,~,~,~,~,nam] = get_obs(dat); % Image (obs) and image properties   

modality = dat.modality{1}.name;             % Imaging modality
I        = size(obs,1);                      % Total number of voxels
ff       = get_ff(vs_s);                     % Fudge factor (see spm_preproc8)   
do_reg   = opt.reg.do_aff || opt.reg.do_nl ; % Update registration?

miss = get_par('missing_struct',obs);

% Bias field
do_bf = opt.bf.do && strcmpi(modality,'MRI'); % Update bias-field?
if do_bf
    % Compute bias-field
    bf = get_bf(dat.bf.chan,dm_s);    
         
    if opt.bf.mc_bf && it_mod > 1
        % Bias field has been mean corrected -> adjust lblnDetbf accordingly        
        lblnDetbf                   = bf;
        lblnDetbf(isnan(lblnDetbf)) = 1;        
        lblnDetbf                   = log(prod(lblnDetbf,2));  
        lblnDetbf                   = sum(lblnDetbf);
        
        dat.lb.lnDetbf(end + 1) = lblnDetbf;
    end
else  
    % Bias field not estimated
    bf = 1;                                    
end

if it_mod > 1
    % Adjust lower bound when building template
    %----------------------------------------------------------------------
    
    append2lb = false;
    
    if opt.reg.mc_aff ...
        % Affine parameters have been mean corrected -> adjust lb accordingly
        aff_reg                 = double(-0.5*dat.reg.r(:)'*opt.reg.aff_reg_ICO*dat.reg.r(:));    
        dat.lb.aff_reg(end + 1) = aff_reg;   

        append2lb = true;
    end

    if opt.model.PropPrior.do && it_mod > opt.reg.strt_nl
        % Tissue proportion prior has been updated
        dat.lb.prop_reg(end + 1) = dot(PropPrior.alpha - 1,log(spm_matcomp('softmax',dat.gmm.prop)));
        
        append2lb = true;
    end
    
    if append2lb
        % Append lower bound
        dat.lb = check_convergence('ini',dat.lb,1,opt.verbose.gmm);
    end
end

% Get initial velocities
%--------------------------------------------------------------------------
if isnumeric(dat.reg.v), v = dat.reg.v;                            
else,                    v = single(dat.reg.v.dat(:,:,:,:));       
end

% FFT of Green's function
prm_v                           = [vs_s ff*opt.reg.rparam];
if opt.reg.int_args > 1, Greens = spm_shoot_greens('kernel',dm_s(1:3),prm_v);
else,                    Greens = [];
end

% Deformation
y = make_deformation(v,prm_v,opt.reg.int_args,Greens);

if opt.verbose.gmm >= 3
    % Verbose stuff
    [~,~,~,~,~,~,~,chn_names] = obs_info(dat);    
    chn_names{end + 1}        = 'Template';
    chn_names{end + 1}        = 'Z';
end

if ~opt.template.do
    % Set initial posteriors to same values as GaussPrior
    dat.gmm.cluster{1}{1} = GaussPrior{1};
    dat.gmm.cluster{1}{2} = GaussPrior{2};
    dat.gmm.cluster{2}{1} = GaussPrior{3};
    dat.gmm.cluster{2}{2} = GaussPrior{4};
end

%--------------------------------------------------------------------------
% Initial alignment of template
%--------------------------------------------------------------------------

if ~opt.template.do && dm_s(3) > 1
    % Align template by mutual information registration
    
    dat = maffreg_template2subject(dat,model,opt);
    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = model.template.nii.mat\E*mat_s;

    % Warp template to subject      
    Template = warp_template(model,y,Affine);
elseif it_mod == 1 && opt.reg.do_aff        
    % Align the template by updating the affine parameters until convergence    

    % Affine matrix    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = mat_a\E*mat_s;

    % Warp template to subject
    Template = warp_template(model,y,Affine);  

    % Disable MRF just here
    do_mrf0    = dat.mrf.do;
    dat.mrf.do = false;
    
    for it_reg=1:opt.reg.nit_init_aff
        
        if opt.verbose.reg >= 3 
            % Show observed data, warped template and current responsibilities
            Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
            show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality,chn_names,opt.model.nam_cls);
            clear Z
        end
        
        [dat,Affine,Template,gain] = update_affine(dat,model,obs,Template,bf,labels,mat_a,mat_s,y,dm_s,scl,miss,opt);
        
        if gain < opt.reg.init_aff_tol
            % Finished updating registration
            break;
        end
    end
    
    dat.mrf.do = do_mrf0;
else
    % Template already aligned
    
    % Affine matrix    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = mat_a\E*mat_s;

    % Warp template to subject
    Template = warp_template(model,y,Affine);  
end

if ~strcmpi(modality,'ct') && ~opt.template.do
    % Get posteriors using the aligned template
    dat.gmm.cluster = get_cluster(obs,bf,dm_s,GaussPrior,miss,{Template,dat.gmm.prop,labels,dat.gmm.part},'sort_pars',false);
end

if opt.do.mrf
    % When using an MRF, store previous responsibilities, etc., in the mrf
    % field of dat
    dat.mrf.oZ = reshape(get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt,'use_uint',true),[dm_s K]); 
    % figure; imshow3D(squeeze(dat.mrf.oZ))
    
    if opt.verbose.mrf >= 3
        show_mrf(dat.mrf);
    end
end

%--------------------------------------------------------------------------
% Verbose
%--------------------------------------------------------------------------
    
if opt.verbose.gmm >= 3
    plot_subject_lowerbound(dat.lb);
end
if opt.verbose.gmm >= 3
    % Show observed data, warped template and current responsibilities
    Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
    show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality,chn_names,opt.model.nam_cls);
    clear Z
end
if opt.verbose.gmm >= 4
    % Show GMM
    show_gmm(dat,obs)
end
if opt.verbose.bf >= 3
    % Show bias field
    show_bf_and_ivel(obs,dm_s,bf);
end
if opt.verbose.reg >= 3 
     % Show initial velocities
    show_bf_and_ivel(obs,dm_s,v);
end
if opt.verbose.reg >= 3 || opt.verbose.gmm >= 3  || opt.verbose.bf >= 3
    % Distribute figures
    deal_figs(model);
end

%--------------------------------------------------------------------------
% Start EM segmentation algorithm
%--------------------------------------------------------------------------

if opt.verbose.gmm, tic; end % Start timer

for it_seg=1:opt.seg.niter     
    
    do_prop = opt.prop.do   && (it_mod >= opt.start_it.do_prop || it_seg >= opt.start_it.do_prop);    
    do_mg   = opt.do.mg     && (it_mod >= opt.start_it.do_mg   || it_seg >= opt.start_it.do_mg);
    do_nl   = opt.reg.do_nl && (it_mod >= opt.reg.strt_nl      || it_seg >= opt.reg.strt_nl);
    
    if it_seg == 1 && do_mg && ~opt.template.do && ~isequal(GaussPrior{7},get_par('lkp',modality,opt))
        % Introduce multiple Gaussians per tissue
        [dat,GaussPrior] = introduce_lkp(dat,model,opt,GaussPrior);                        
    end
    
    %----------------------------------------------------------------------
    % UPDATE BIAS FIELD
    %----------------------------------------------------------------------
            
    for it_bf=1:opt.bf.niter                 
        % Start with updating GMM parameters
        dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,do_mg,opt);                   
                
        if (it_bf > 1 && ~((dat.lb.last - ooll)>2*opt.bf.tol*I)) || it_bf == opt.bf.niter, break; end
        ooll = dat.lb.last;
                
        if do_bf 
            % Update bias-field parameters    
            [dat,bf] = update_bf(dat,obs,bf,Template,dm_s,labels,miss,opt);             
        end               
    
        if opt.verbose.gmm >= 3                        
            show_seg(obs,Template,dat.gmm.prop, ...
                     get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                     dm_s,modality,chn_names,opt.model.nam_cls);            
        end                
    end   
       
    %----------------------------------------------------------------------
    % UPDATE REGISTRATION
    %----------------------------------------------------------------------
        
    if do_reg              
        for it_reg=1:opt.reg.niter 
        
            if opt.reg.do_aff                
                % Update affine parameters                                
                [dat,Affine,Template,gain] = update_affine(dat,model,obs,Template,bf,labels,mat_a,mat_s,y,dm_s,scl,miss,opt);               
            
                if opt.verbose.gmm >= 3
                    show_seg(obs,Template,dat.gmm.prop, ...
                             get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                             dm_s,modality,chn_names,opt.model.nam_cls);        
                end   
            end
            
            if do_nl
                % Update initial velocities                                                  
                [dat,y,v,Template,Greens,gain,opt] = update_nonlin(dat,model,obs,Template,bf,labels,v,y,Affine,Greens,it_seg,it_mod,scl,miss,opt);                

                if opt.verbose.gmm >= 3
                    show_seg(obs,Template,dat.gmm.prop, ...
                             get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                             dm_s,modality,chn_names,opt.model.nam_cls);        
                end    
            end        
            
            if gain < opt.reg.tol
               % Finished updating registration
               break;
            end
        end
    end  
    
    %----------------------------------------------------------------------
    % UPDATE TISSUE PROPORTIONS
    %----------------------------------------------------------------------
    
    if do_prop         
        for it_bf=1:opt.bf.niter                 
            % Start with updating GMM parameters
            dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,do_mg,opt);                   

            if (it_bf > 1 && ~((dat.lb.last - oll)>2*opt.bf.tol*I)) || it_bf == opt.bf.niter, break; end
            oll = dat.lb.last;

            % Update tissue proportions
            dat = update_prop(obs,bf,dat,Template,labels,scl,dm_s,PropPrior,miss,opt);          

            if opt.verbose.gmm >= 3                        
                show_seg(obs,Template,dat.gmm.prop, ...
                         get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                         dm_s,modality,chn_names,opt.model.nam_cls);            
            end
        end   
    end
    
    %----------------------------------------------------------------------
    % UPDATE MRF WEIGHTS
    %----------------------------------------------------------------------
    
    if opt.do.update_mrf && opt.do.mrf && (it_mod >= opt.start_it.do_upd_mrf || it_seg >= opt.start_it.do_upd_mrf)
        for it_bf=1:2
            % Start with updating GMM parameters
            dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,do_mg,opt);                   

            if it_bf == 2, break; end        

            % Update MRF weights        
            dat = update_mrf(dat,obs,bf,Template,dm_s,labels,scl,miss,opt);                    

            if opt.verbose.gmm >= 3                        
                show_seg(obs,Template,dat.gmm.prop, ...
                         get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                         dm_s,modality,chn_names,opt.model.nam_cls);            
            end
        end  
    end
    
    %----------------------------------------------------------------------
    % CONVERGED?
    %----------------------------------------------------------------------
    
    if opt.verbose.gmm                
        % Print overall convergence        
        fprintf('%3s | %3d | lb = %10.6f | diff = %1.5f | tol = %1.5f\n', 'seg', it_seg, dat.lb.last, gain, 1e-4);
    end
    
    if it_seg >= 10 && ~((dat.lb.last - ooll) > 2*opt.seg.tol*I)
        if opt.verbose.gmm
            fprintf('Segmentation converged in %i iterations and %0.0f seconds.\n',it_seg,toc);
        end
        
        break
    end
end
clear Greens

if opt.verbose.model >= 3
   % Write 2D versions of responsibilities to disk (for verbose) 
   ix_z = floor(dm_s(3)/2) + 1;
   Z    = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt,'z',ix_z);
   
   dat.pth.seg2d = fullfile(opt.dir_seg2d,['seg2d_' nam '.nii']);
   spm_misc('create_nii',dat.pth.seg2d,Z,mat_s,[spm_type('float32') spm_platform('bigend')],'seg2d');        
   clear Z
end

if do_nl
    % Save initial velocities
    if isnumeric(dat.reg.v), dat.reg.v              = v;
    else,                    dat.reg.v.dat(:,:,:,:) = v;
    end
end
clear v

if    opt.template.do && opt.template.load_a_der
    %----------------------------------------------------------------------
    % Compute template derivatives (will be loaded later in
    % update_template)
    %----------------------------------------------------------------------

    % Get responsibilities
    Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);   
    clear Template

    % Push responsibilities from subject space to template space
    y = spm_warps('transform',Affine,y);
    if dm_s(3) == 1
        y(:,:,:,3) = 1;
    end
    
    Z = push_responsibilities(Z,y,dm_a(1:3));                      
    clear y
    
    % Compute template gradients and Hessian
    a               = rotate_template(model.template.nii,opt);
    [gr,H,dll]      = diff_template(a,Z,dat.gmm.prop,opt);         
    dat.template.ll = dll;                        
    clear a Z
    
    % Write derivatives to disk                                    
    spm_misc('create_nii',dat.template.pth_gr,gr,mat_a,[spm_type('float32') spm_platform('bigend')],'a_gr'); 
    spm_misc('create_nii',dat.template.pth_H,H,mat_a,[spm_type('float32') spm_platform('bigend')],'a_H');  
    clear gr H
end

if isfield(dat.mrf,'oZ') && opt.template.do
    % Remove oZ field from dat.mrf
    dat.mrf = rmfield(dat.mrf,'oZ');
end

if ~opt.template.do || (opt.template.do && it_mod == opt.model.niter)
    %----------------------------------------------------------------------
    % Write resulting segmentations (etc) to disk
    %----------------------------------------------------------------------
    
    res = write_results(dat,model,opt);
else
    res = [];
end
%==========================================================================