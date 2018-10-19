function Z = clean_eyeballs(Z,model,y,opt)

% Get template (in default space)
Nii       = model.template.nii;
Template  = single(Nii.dat(:,:,:,:));
Template  = spm_matcomp('softmax',Template);  
K         = model.template.nii.dat.dim(4);

% Get background and foreground classes
[bg0,fg0] = get_par('bg_fg',opt.seg.bg,K);
    
% Get background mask from template
bg  = sum(Template(:,:,:,bg0),4);
msk = bg >= 0.1; 
clear bg Template

% Grow background a bit into brain region by dilation
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;
th    = 0.6;
niter = 32;
for j=1:niter        
    b = single(msk > th);

    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end
msk = b > 0;

% 'Inverse' bacground mask to get brain mask
msk = ~msk;

% Find largest connected component and fill holes
L   = bwlabeln(msk);
msk = L == 1;
msk = imfill(msk,'holes');

% Warp brain mask to subject
msk                 = spm_diffeo('bsplins',single(msk),y,[0 0 0  0 0 0]);
msk(~isfinite(msk)) = 0;

% Clean-up brain classes of responsibilities with constructed brain mask
Z(:,:,:,fg0) = bsxfun(@times,Z(:,:,:,fg0),msk);
Z            = bsxfun(@rdivide,Z,sum(Z,4) + eps);
%==========================================================================