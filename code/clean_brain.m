function [Z,msk_les] = clean_brain(Z,model,y,opt)
% FORMAT [Z,msk_les] = clean_brain(Z,model,y,opt)
% Z       - Class responsibilities [Nx Ny Nz K]
% model   - Model structure
% y       - Transform to warp template to subject [Nx Ny Nz 3]
% opt     - Options structure
% msk_les - Lesion mask
%
% Clean-up resulting (soft) segmentations using a series of educated
% heuristics.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
speak = false;

P0 = single(model.template.nii.dat(:,:,:,:));
P0 = spm_matcomp('softmax',P0);  

% Get background and foreground classes
K         = size(P0,4);
[bg0,fg0] = get_par('bg_fg',opt.seg.bg,K);
    
% Get background mask from template
P = sum(P0(:,:,:,bg0),4);

if speak
    dm        = size(P0);
    dm0       = size(Z);
    step      = ceil(0.05*dm(3));
    step0     = ceil(0.05*dm0(3));

    figure(101);
    subplot(131)
    img = [];
    for i=1:step:dm(3)
        img = [img; P(:,:,i)];
    end
    imagesc(img); axis off xy
    drawnow
end

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Dilate background mask
% dilate (bright zones gets bigger) if <0.5, erode (dark zones gets bigger) if th>0.5
%--------------------------------------------------------------------------
b      = P;
th     = 0.4; 
thrsld = 0.1;
niter  = 128;
for j=1:niter
    
    for i=1:size(b,3)
        bg              = double(P(:,:,i));
        bg(bg > thrsld) = 1;
        
        bp       = double(b(:,:,i));
        bp       = (bp>th).*(bg);
        b(:,:,i) = bp;
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);  
    
    if speak
        subplot(132)
        img = [];
        for i=1:step:dm(3)
            img = [img; b(:,:,i)];
        end
        imagesc(img); axis off xy
        drawnow
    end
end

% Get brain mask
%--------------------------------------------------------------------------
msk = b > thrsld;
msk = ~msk;

% Find largest connected component in brain mask
%--------------------------------------------------------------------------
L   = bwlabeln(msk);
nL  = unique(L);
vls = zeros([1 numel(nL)]);
for i=1:numel(nL)
    vls(i) = sum(sum(sum(L == nL(i))));
end
[~,ix] = sort(vls);
ix     = (ix(end - 1) - 1);
msk    = L == ix;

if speak
    subplot(133)
    img = [];
    for i=1:step:dm(3)
        img = [img; msk(:,:,i)];
    end
    imagesc(img); axis off xy
    drawnow
end

% Warp brain mask to subject
%--------------------------------------------------------------------------
msk                 = spm_diffeo('bsplins',single(msk),y,[0 0 0  0 0 0]);
msk(~isfinite(msk)) = 0;

% Clean-up brain classes of responsibilities with constructed brain mask
%--------------------------------------------------------------------------
Z(:,:,:,fg0) = bsxfun(@times,Z(:,:,:,fg0),msk);
Z            = bsxfun(@rdivide,Z,sum(Z,4) + eps);

if speak
    figure(102)
    for k=1:K
        subplot(1,K,k)
        img = [];
        for i=1:step0:dm0(3)
            R = Z(:,:,i,k);    
            if find(k == fg0)
                R(~msk(:,:,i)) = 0;
            end
            img = [img; R];
        end
        imagesc(img); axis off xy
        drawnow
    end
end

% Erode brain mask to create lesion mask
%--------------------------------------------------------------------------
b     = single(msk);
th    = 0.6;
niter = 20;
for j=1:niter
    
    for i=1:size(b,3)
        bg              = double(msk(:,:,i));
        
        bp       = double(b(:,:,i));
        bp       = (bp>th).*(bg);
        b(:,:,i) = bp;
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);  
    
    if speak
        figure(101)
        subplot(133)
        img = [];
        for i=1:step0:dm0(3)
            R = Z(:,:,i,4);  
            msk0 = b(:,:,i) > 0;
            R(~msk0) = 0;
            img = [img; R];
        end
        imagesc(img); axis off xy
        drawnow
    end
end

msk_les = b > 0; % Lesion mask
%==========================================================================