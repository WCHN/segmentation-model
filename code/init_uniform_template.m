function model = init_uniform_template(dat,opt)
K         = opt.template.K;
dir_model = opt.dir_model;
S0        = numel(dat);
vs        = opt.template.vs;

% Collect orientation matrices and dimensions from all input images
%--------------------------------------------------------------------------
mats = [];
dms  = [];   
for s=1:S0
    if isfield(dat{s}.modality{1},'channel')
        C = numel(dat{s}.modality{1}.channel);
        for c=1:C
            mat  = dat{s}.modality{1}.channel{c}.nii.mat;           
            mats = cat(3,mats,mat);

            dm  = dat{s}.modality{1}.channel{c}.nii.dat.dim;                          
            if numel(dm)==2, dm(3) = 1; end
            dms = cat(1,dms,dm);
        end
    else
        mat  = dat{s}.modality{1}.nii.mat;           
        mats = cat(3,mats,mat);

        dm  = dat{s}.modality{1}.nii.dat.dim;                          
        if numel(dm)==2, dm(3) = 1; end
        dms = cat(1,dms,dm);
    end
end

% Compute average orientation matrix and dimensions
[M,d] = spm_misc('compute_avg_mat',mats,dms);

% Write template to disk
%--------------------------------------------------------------------------
pth_template  = fullfile(dir_model,'template.nii');   
[pth,nam,ext] = fileparts(pth_template);

vols = cell(K,1);
for k=1:K    
    
    img = zeros(d,'single');
    
    vols{k} = fullfile(pth,[nam num2str(k) ext]);
    spm_misc('create_nii',vols{k},img,M,[spm_type('float32') spm_platform('bigend')],'template');
end
clear img

matlabbatch                       = cell(1,1);
matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_template;    
matlabbatch{1}.spm.util.cat.dtype = 0;
spm_jobman('run',matlabbatch);

for k=1:K
    delete(vols{k});
end

% Set template
model              = struct;
model.template.nii = nifti(pth_template);

% Crop template and/or change voxel sizes
model = resize_template(model,opt);
d     = model.template.nii.dat.dim;

for k=1:K        
    img = log(ones(d(1:3),'single')/K);
    
    if ~isempty(opt.lesion.hemi)
        % Assumes there are two lesion classes, one for each hemisphere
        if k == opt.lesion.hemi{1}
            img(1:floor(d(1)/2 - 0.05*d(1)),:,:)  = log(1e-2);
            img = spm_imbasics('smooth_img',img,12);
        elseif k == opt.lesion.hemi{2}
            img(ceil(d(1)/2 + 0.05*d(1)):end,:,:) = log(1e-2);
            img = spm_imbasics('smooth_img',img,12);
        end
    end
    
    model.template.nii.dat(:,:,:,k) = img;
end

% Get values to be used for FOV voxels when warping 
model = init_template_bg(model,opt);

% Init template objval
model.template.objval.post  = 0;
model.template.objval.likel = 0;
model.template.objval.pr    = 0;

alpha           = ones(1,K)*opt.prop.reg;
PropPrior.alpha = alpha;
PropPrior.norm  = 0;
model.PropPrior = PropPrior;

% Save PropPrior
fname = fullfile(opt.dir_model,'PropPrior.mat');
save(fname,'PropPrior');
%==========================================================================