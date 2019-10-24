function [P] = clean_gwc(P,ix,level)
if nargin < 2
    % Default SPM12 template ordering
    ix.gm = 1;
    ix.wm = 2;
    ix.cs = 3;
end
if nargin < 3, level = 1; end

b = sum(P(:,:,:,ix.wm),4);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(sum(P(:,:,i,ix.gm),4));
        wp       = double(sum(P(:,:,i,ix.wm),4));
        bp       = double(b(:,:,i));
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = bp;
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end

% Also clean up the CSF.
if niter2 > 0
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = double(sum(P(:,:,i,ix.gm),4));
            wp       = double(sum(P(:,:,i,ix.wm),4));
            cp       = double(sum(P(:,:,i,ix.cs),4));
            bp       = double(c(:,:,i));
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = bp;
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4)
        slices{k1} = double(P(:,:,i,k1));
    end
    bp           = double(b(:,:,i));
    bp           = ((bp>th).*(sum(cat(3,slices{ix.gm}),3)+sum(cat(3,slices{ix.wm}),3)))>th;
    for i1=1:numel(ix.gm)
        slices{ix.gm(i1)} = slices{ix.gm(i1)}.*bp;
    end
    for i1=1:numel(ix.wm)
        slices{ix.wm(i1)} = slices{ix.wm(i1)}.*bp;
    end
    
    if niter2>0
        cp           = double(c(:,:,i));
        cp           = ((cp>th).*(sum(cat(3,slices{ix.gm}),3)+sum(cat(3,slices{ix.wm}),3)+sum(cat(3,slices{ix.cs}),3)))>th;
        
        for i1=1:numel(ix.cs)
            slices{ix.cs(i1)} = slices{ix.cs(i1)}.*cp;
        end        
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4)
        P(:,:,i,k1) = slices{k1}./tot;
    end 
end