function res = write_results(dat,model,opt)

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

[obs,dm_s,mat_s,vs_s,scl]  = get_obs(dat); 
[~,~,~,C,nam,~,obs_fnames] = obs_info(dat);
miss                       = get_par('missing_struct',obs);
dm_a                       = model.template.nii.dat.dim;
mat_a                      = model.template.nii.mat;
vs_a                       = sqrt(sum(mat_a(1:3,1:3).^2));
K                          = dm_a(4);
modality                   = dat.modality{1}.name; 
ff                         = get_ff(vs_s);   
if isnumeric(dat.reg.v)
    v                      = dat.reg.v;                       
else
    v                      = single(dat.reg.v.dat(:,:,:,:));   
end
labels                     = get_labels(dat,opt);
prm_v                      = [vs_s ff*opt.reg.rparam];   
do_bf                      = opt.bf.do && strcmpi(modality,'MRI');
if do_bf, bf               = get_bf(dat.bf.chan,dm_s);
else,     bf               = 1;
end

%--------------------------------------------------------------------------
% Prepare function output
%--------------------------------------------------------------------------
res    = struct;
res.c  = cell(K,1);
res.wc = cell(K,1);

%--------------------------------------------------------------------------
% Get deformation
%--------------------------------------------------------------------------

if opt.reg.int_args > 1, Greens = spm_shoot_greens('kernel',dm_s(1:3),prm_v);
else,                    Greens = [];
end

y = make_deformation(v,prm_v,opt.reg.int_args,Greens);  

%--------------------------------------------------------------------------
% Write initial velocity (v)
%--------------------------------------------------------------------------

Nii         = nifti;
Nii.dat     = file_array(fullfile(dat.dir.def,['vel_', nam{1}, '.nii']),...
                         size(v),...
                         [spm_type('float32') spm_platform('bigend')],...
                         0,1,0);
Nii.mat     = mat_s;
Nii.mat0    = mat_s;
Nii.descrip = 'Initial velocity';
create(Nii);

Nii.dat(:,:,:,:) = v;

%--------------------------------------------------------------------------
% Warp template to subject
%--------------------------------------------------------------------------

E        = spm_dexpm(dat.reg.r,opt.reg.B); 
Affine   = model.template.nii.mat\E*mat_s;               
Template = warp_template(model,y,Affine);
y        = spm_warps('transform',Affine,y);
if dm_s(3) == 1
    y(:,:,:,3) = 1;
end

%--------------------------------------------------------------------------
% Get responsibilities
%--------------------------------------------------------------------------

Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
clear Template labels

%--------------------------------------------------------------------------
% Write non-preprocessed responsibilities to disk
%--------------------------------------------------------------------------

for k=1:size(Z,4)
    fname = fullfile(dat.dir.seg_orig,['resp' num2str(k) nam{1} '.nii']);
    spm_misc('create_nii',fname,Z(:,:,:,k),mat_s,[spm_type('float32') spm_platform('bigend')],'Segmentation');
end

%--------------------------------------------------------------------------
% Write images (bias-field corrected)
%--------------------------------------------------------------------------

for c=1:C    
    Nii         = nifti;
    Nii.dat     = file_array(fullfile(dat.dir.img,['Corrected_', nam{c}, '.nii']),...
                             dm_s,...
                             [spm_type('float32') spm_platform('bigend')],...
                             0,1,0);
    Nii.mat     = mat_s;
    Nii.mat0    = mat_s;
    Nii.descrip = 'Bias Field Corrected Image';
    create(Nii);
    
    Nii.dat(:,:,:) = reshape(obs(:,c).*bf(:,c),dm_s);
end

%--------------------------------------------------------------------------
% Write bias field
%--------------------------------------------------------------------------
    
if opt.bf.do && strcmpi(modality,'MRI')
    for c=1:C
        Nii         = nifti;
        Nii.dat     = file_array(fullfile(dat.dir.bf,['BiasField_', nam{c}, '.nii']),...
                                 dm_s,...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        Nii.mat     = mat_s;
        Nii.mat0    = mat_s;
        Nii.descrip = 'Estimated Bias Field';
        create(Nii);

        Nii.dat(:,:,:) = reshape(bf(:,c),dm_s);
    end
end

%--------------------------------------------------------------------------
% Compute inverse deformation (iy)
%--------------------------------------------------------------------------

[~,~,~,iy] = spm_shoot3d(v,prm_v,[opt.reg.int_args [2 2]],Greens);
iy         = spm_warps('compose',iy,inv(Affine),single(spm_warps('identity',dm_a(1:3))));
clear Greens v

if dm_s(3) == 1
    iy(:,:,:,3) = 1;
end

%--------------------------------------------------------------------------
% Write bias-field corrected images, warped to MNI space
%--------------------------------------------------------------------------

for c=1:C    
    Nii         = nifti;
    Nii.dat     = file_array(fullfile(dat.dir.img,['MNI_Corrected_', nam{c}, '.nii']),...
                             dm_a(1:3),...
                             [spm_type('float32') spm_platform('bigend')],...
                             0,1,0);
    Nii.mat     = mat_a;
    Nii.mat0    = mat_a;
    Nii.descrip = 'Bias Field Corrected Image in MNI space';
    create(Nii);
    
    x               = spm_diffeo('bsplins', reshape(obs(:,c).*bf(:,c),dm_s), iy, [1 1 1 0 0 0]);
    x(~isfinite(x)) = 0;
    
    Nii.dat(:,:,:) = single(x);
end
clear bf obs x iy
    
%--------------------------------------------------------------------------
% Some cleaning up of responsibilities
%--------------------------------------------------------------------------

if opt.clean.eyeballs
    % Ad-hoc clean-up to remove crap outside of brain from brain classes           
    Z = clean_eyeballs(Z,model,y,opt);    
    % k=4;d=ceil(dm_s(3)/2);figure(666);subplot(121);imagesc(squeeze(Z(:,:,d,k)));subplot(122);imagesc(squeeze(Z1(:,:,d,k)))
end

if opt.clean.cnn_mrf.do
    % Ad-hoc CNN-MRF clean-up of lesion class (and writing to disk)        
    x    = cell(1,2);
    x{1} = single(Z(:,:,:,opt.clean.cnn_mrf.les_cl) > opt.clean.cnn_mrf.val);
    x{2} = 1 - x{1};

    im_clean = spm_cnn_mrf('predict',opt.clean.cnn_mrf.pth_net,x,'speak',0,'w_on',0.9,'nit',40);
    im_clean = im_clean{1};
    clear x
    
    % Get largest connected component
    les_msk = im_clean > opt.clean.cnn_mrf.val;
    L       = bwlabeln(les_msk);
    nL      = unique(L);
    vls     = zeros([1 numel(nL)]);
    for i=1:numel(nL)
        vls(i) = sum(sum(sum(L == nL(i))));
    end
    [~,ix] = sort(vls);
    
    les_msk = L == (ix(end - 1) - 1);
    les_msk = imfill(les_msk,'holes');
%     figure;imshow3D(les_msk)    
    les_msk = single(les_msk);
    clear im_clean
    
    fname = fullfile(dat.dir.seg,['cles_' nam{1} '.nii']);
    spm_misc('create_nii',fname,les_msk,mat_s,[spm_type('float32') spm_platform('bigend')],'Lesion (native)');
        
    res.cles = nifti(fname);
    
    if 0 && isfield(dat,'label')
        gt  = dat.label{1}.nii.dat(:,:,:) > 0;
        prd = logical(les_msk);        
        fprintf('ds=%0.3f\n',dice(prd,gt));
    end
end

if opt.clean.mrf.do
    % Ad-hoc MRF clean-up of segmentation
    Z = mrf_post(Z,vs_s,opt);
end

%--------------------------------------------------------------------------
% Write final segmentation to disk
%--------------------------------------------------------------------------

% Native space (c)
for k=1:size(Z,4)
    fname = fullfile(dat.dir.seg,['c' num2str(k) nam{1} '.nii']);
    spm_misc('create_nii',fname,Z(:,:,:,k),mat_s,[spm_type('float32') spm_platform('bigend')],'Segmentation (native)');
    
    res.c{k} = nifti(fname);
end

% Warped to template space (wc)   
C = zeros(dm_a,'single');
for k=1:size(Z,4)
    [c,w]      = spm_diffeo('push',Z(:,:,:,k),y,dm_a(1:3));    
    C(:,:,:,k) = spm_field(w,c,[vs_a  1e-6 1e-4 0  3 2]);
    clear w c
end

C = max(C,eps);
s = sum(C,4);
for k=1:size(Z,4)    
    fname = fullfile(dat.dir.seg,['wc' num2str(k) nam{1} '.nii']);
    spm_misc('create_nii',fname,C(:,:,:,k)./s,mat_a,[spm_type('float32') spm_platform('bigend')],'Segmentation (template)');
    
    res.wc{k} = nifti(fname);
end
clear s C

if opt.clean.cnn_mrf.do
    [c,w]   = spm_diffeo('push',les_msk,y,dm_a(1:3));    
    les_msk = spm_field(w,c,[vs_a  1e-6 1e-4 0  3 2]);
    clear w c
    
    fname = fullfile(dat.dir.seg,['wcles_' nam{1} '.nii']);
    spm_misc('create_nii',fname,les_msk > opt.clean.cnn_mrf.val,mat_a,[spm_type('float32') spm_platform('bigend')],'Lesion (template)');
    clear les_msk
    
    res.wcles = nifti(fname);
end
clear y

%--------------------------------------------------------------------------    
% Write ML labels to disk
%--------------------------------------------------------------------------
    
if opt.seg.write_mllabels    
    [bg,fg] = get_par('bg_fg',opt.seg.bg,K);    
    Z       = cat(4,sum(Z(:,:,:,bg),4),Z(:,:,:,fg));    

    [~,ml] = max(Z,[],4);
    ml     = ml - 1;

    fname = fullfile(dat.dir.ml,'ml.nii');
    spm_misc('create_nii',fname,ml,mat_s,[spm_type('float32') spm_platform('bigend')],'ML-labels');    
end
        
%--------------------------------------------------------------------------    
% Show final segmentations
%--------------------------------------------------------------------------

if opt.seg.show
    files = spm_select('FPList',dat.dir.seg,'^c.*\.nii$');
    spm_check_registration(char(obs_fnames, files));
end

%--------------------------------------------------------------------------    
% Clean-up temporary files
%--------------------------------------------------------------------------

if opt.model.clean_up && isfield(dat,'dir') && isfield(dat.dir,'vel')
    rmdir(dat.dir.vel,'s');
end
%==========================================================================