function Z = clean_eyeballs(Z,model,y)
Nii      = model.template.nii;
Template = single(Nii.dat(:,:,:,:));
Template = spm_matcomp('softmax',Template);  
tiny     = 1e-4;
K        = model.template.nii.dat.dim(4);

[bg,fg] = get_par('bg_fg',K);
    
bg  = sum(Template(:,:,:,bg),4);
clear Template

msk = bg >= 0.1; 
clear bg

msk = ~msk;

L   = bwlabeln(msk);
msk = L == 1;
msk = imfill(msk,'holes');

% Warp mask to subject
msk                 = spm_diffeo('bsplins',single(msk),y,[0 0 0  0 0 0]);
msk(~isfinite(msk)) = 0;

if 0
%     figure(667); imshow3D(msk);    
end

% % Dilate 
% kx=[0.75 1 0.75];
% ky=[0.75 1 0.75];
% kz=[0.75 1 0.75];
% sm=sum(kron(kron(kz,ky),kx))^(1/3);
% kx=kx/sm; ky=ky/sm; kz=kz/sm;
% th   = 0.6;
% niter = 32;
% for j=1:niter        
%     b = single(msk > th);
% 
%     spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
% end
% msk = b > 0;

if 0
%     figure(667); imshow3D(msk);    
end

Z(:,:,:,fg) = bsxfun(@times,Z(:,:,:,fg),msk);
Z           = bsxfun(@rdivide,Z,sum(Z,4) + eps);

if 0    
    nz = floor(size(Z,3)/2) + 1;
%     figure(667); imshow3D(squeeze(Z(:,:,nz,:)));
end
%==========================================================================