function [model,dat] = update_template(dat,model,opt,is_init)
% FORMAT [model,dat] = update_template(dat,model,opt,is_init)
% dat     - Subject data structure
% model   - Model structure
% opt     - Options structure
% is_init - 
%
% Update the template by Gauss-Newton.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4, is_init = false; end
if is_init
    opt.template.niter = 1;             % Number of GN iterations
end

% Options
reg        = opt.template.reg;          % [a m b] Template regularisation
verbose    = opt.template.verbose;      % Verbosity level
shrink     = opt.template.shrink;       % Crop template so it fits pushed data
B          = opt.reg.B;                 % Affine Lie basis
prm_reg    = opt.reg.rparam;            % [a m b le1 le1] Velocity regularisation
int_args   = opt.reg.int_args;          % Number of integration steps
load_a_der = opt.template.load_a_der;   % Precomputed template derivatives
R          = opt.template.R;            % Template null-space rotation
nii_a      = model.template.nii;        % Nifti holding the template

% Parameters
S0     = numel(dat);                    % Number of subjects
dm     = [nii_a.dat.dim,1,1,1];        
K      = dm(4);                         % Number of template classes
rits   = [3 3];                         % FMG cycles and relaxation iterations
mat_a  = nii_a.mat;                     % Template vox2world matrix
dm_a   = nii_a.dat.dim;                 % Template dimensions
vs_a   = spm_misc('vxsize',mat_a);      % Template voxel size
prm_a  = [vs_a prod(vs_a)*reg];         % Template reg param vector

% Initial starting estimates
a = rotate_template(nii_a,opt);         % Send previous template in null space

%--------------------------------------------------------------------------
% Compute/load pushed responsibilities then update template using Gauss-Newton
%--------------------------------------------------------------------------

if load_a_der && ~is_init
    % Template derivatives are loaded from disk
    %----------------------------------------------------------------------
                    
    % Add to global derivatives using bounding boxes
    H  = zeros([dm(1:3) round(((K-1)*K)/2)],'single'); % 2nd derivatives
    gr = zeros([dm(1:3),K-1],'single');                % 1st derivatives
    ll = 0;
    for s=1:S0
        % Add up gradients
        nii = nifti(dat{s}.template.pth_gr);
        gr  = gr + nii.dat(:,:,:,:);
        
        % Add up Hessians
        nii = nifti(dat{s}.template.pth_H);
        H   = H + nii.dat(:,:,:,:);
        
        % Add up ll:s
        ll = ll + dat{s}.template.ll;
    end            
    
    % Do Gauss-Newton step to update template
    [a,~,ll1,ss2] = solve_template(a,gr,H,prm_a,rits);
    
    if verbose>=2
        fprintf('update_template | %2d | ll = %6.6f    ll_a = %6.6f    ll + ll_a = %6.6f    ss2 = %6.6f\n',1,ll/prod(dm(1:3)), ll1/prod(dm(1:3)), (ll+ll1)/prod(dm(1:3)), (ss2)/prod(dm(1:3)));
    end
else
    % Template derivatives are computed on the fly.
    % This allows for iterating the Gauss-Newton solver.
    %----------------------------------------------------------------------
    
    for iter=1:opt.template.niter
        H  = zeros([dm(1:3) round(((K-1)*K)/2)],'single'); % 2nd derivatives
        gr = zeros([dm(1:3),K-1],'single');                % 1st derivatives    
        ll = 0;
%         for s=1:S0
        parfor s=1:S0 
            
            % Subject parameters
            [obs,dm_s,mat_s,vs_s,scl,~,~,~,~,nam] = get_obs(dat{s});                
            
            ff                        = get_ff(vs_s);     
            prop                      = dat{s}.gmm.prop;
            prm_v                     = [vs_s ff*prm_reg];
            if isnumeric(dat{s}.reg.v)
                v                     = dat{s}.reg.v;                            % Initial velocity stored in array
            else
                v                     = single(dat{s}.reg.v.dat(:,:,:,:));       % Initial velocity read from NIfTI
            end            
            modality                  = dat{s}.modality{1}.name; 
            do_bf                     = opt.bf.do && strcmpi(modality,'MRI');
            if do_bf || strcmpi(modality,'MRI') 
                bf              = get_bf(dat{s}.bf.chan,dm_s);
            else     
                bf              = 1;
            end
            labels                    = get_labels(dat{s},opt);
            miss                      = get_par('missing_struct',obs);

            % Build and apply FFT of Green's function to map from momentum
            % to velocity (if opt.reg.int_args > 1)                       
            if int_args > 1, Greens = spm_shoot_greens('kernel',dm_s(1:3),prm_v);
            else             Greens = [];
            end

            % Generates deformations from initial velocity fields by
            % gedesic shooting (if opt.reg.int_args > 1)
            y      = make_deformation(v,prm_v,int_args,Greens);
            Greens = [];
            v      = [];
            
            % Compute affine transformation matrix
            E      = spm_dexpm(dat{s}.reg.r,B); % Compute matrix exponential
            Affine = mat_a\E*mat_s;   

            % Warp template to subject and softmax       
            [Template,y] = warp_template(model,y,Affine);              

            % Get responsibilities
            Z        = get_resp(obs,bf,dat{s},Template,labels,scl,miss,dm_s,opt);   
            labels   = []; 
            Template = []; 
            % figure; imshow3D(squeeze(reshape(Z,[dm_s 9])))                                          

            if opt.verbose.model >= 3
                % Write 2D versions to disk (for verbose) of..
                ix_z = floor(dm_s(3)/2) + 1;

                % ..responsibilities 
                dat{s}.pth.seg2d = fullfile(opt.dir_seg2d,['seg2d_' nam '.nii']);
                spm_misc('create_nii',dat{s}.pth.seg2d,Z(:,:,ix_z,:),mat_s,[spm_type('float32') spm_platform('bigend')],'seg2d');                        

                % ..of image (only one channel)
                [~,~,~,~,~,~,~,chn_names] = obs_info(dat{s}); 
                for i=1:numel(chn_names)
                   if strcmpi(chn_names{i},'T1')
                       break
                   end
                end

                im              = reshape(obs(:,i),dm_s(1:3));
                im              = im(:,:,ix_z);
                dat{s}.pth.im2d = fullfile(opt.dir_seg2d,['im2d_' nam '.nii']);
                spm_misc('create_nii',dat{s}.pth.im2d,im,mat_s,[spm_type('float32') spm_platform('bigend')],'im2d');   
                
                dat{s}.pth.bfim2d = fullfile(opt.dir_seg2d,['bfim2d_' nam '.nii']);
                spm_misc('create_nii',dat{s}.pth.bfim2d,im,mat_s,[spm_type('float32') spm_platform('bigend')],'bfim2d');    
                im              = [];
            end            
            
            % Push responsibilities in subject space to template space
            Z = push_responsibilities(Z,y,dm_a(1:3)); 
            y = []; 
            % figure; imshow3D(squeeze(Z))
    
            % Compute gradients and Hessian
            [gr_s,H_s,ll_s] = diff_template(a,Z,prop,opt); 
            Z               = [];            

            % Add to global derivatives using bounding box         
            gr = gr + gr_s;
            H  = H  + H_s;            
            ll = ll + ll_s;                        
        end    

        % Do Gauss-Newton step to update template
        [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm_a,rits);

        if verbose>=2
            fprintf('update_template | %2d | ll = %6.6f    ll_a = %6.6f    ll + ll_a = %6.6f    ss2 = %6.6f\n', iter,ll/prod(dm(1:3)), ll1/prod(dm(1:3)), (ll+ll1)/prod(dm(1:3)), (ss2)/prod(dm(1:3)));
        end

        if ss2/ss1<1e-4 
            % Converged?
            break; 
        end        
    end
end
ll = -(ll + ll1);

if ~isfinite(ll)
    % Just a simple sanity check
    error('~isfinite(ll)');
end

model.template.objval.post(end + 1)  = -(ll + ll1);
model.template.objval.likel(end + 1) = -ll;
model.template.objval.pr(end + 1)    = -ll1;

% Rotate back null-space (soft-maxed version is just used for
% visualisation)
[~,a] = softmax_template(a,R);

% Update the NIfTI that stores the template
nii_a.dat(:,:,:,:) = a;

if shrink
    % Use all the computed bounding-boxes from the push operations to
    % shrink the template.
    nii_a = spm_impreproc('mult_bb_crop',nii_a,BB);
end

model.template.nii = nii_a;

% Update values to be used for FOV voxels when warping 
model = init_template_bg(model,opt);

if verbose>=1
    % Show soft-maxed template
    show_template(model,opt,S0);
end
%==========================================================================

%==========================================================================
function [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm,its)
% FORMAT [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm,its)
% a   - log-template in null space
% gr  - Template gradient
% H   - Template Hessian
% prm - Regularisation parameters  [vs abs memb bend]
% its - Full multi grid parameters [cycles relax] 
%
% Perform one Gauss-Newton update of the template.
% ss1 and ss2 are for examining how close the 1st derivatives are to zero.
% At convergence, the derivatives from the likelihood term should match those
% from the prior (regularisation) term.
K   = size(a,4);
ss1 = sum(sum(sum(sum(gr.^2))));
gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
ll1 = 0;
for k=1:K
    ll1 = ll1 + 0.5*sum(sum(sum(sum(double(gr1(:,:,:,k)).*double(a(:,:,:,k)))))); % -ve log probability of the prior term
end
gr  = gr + gr1;                       % Combine the derivatives of the two terms
clear gr1

if ~isfinite(ll1)
    warning('~isfinite(ll1)');
end
            
ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence

a = a - spm_field(H,gr,[prm(1:3) prm(4) prm(5:6) its]); % Gauss-Newton update  
%==========================================================================