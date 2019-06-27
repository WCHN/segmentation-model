function [P] = clean_gwc(P,level)
if nargin<2, level = 1; end
% GM=1 -> GM=4
% WM=2 -> WM=5
% CS=3 -> CS=3
ixwm = 5;
ixgm = 4;
ixcs = 3;

b    = P(:,:,:,ixwm);

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
        gp       = double(P(:,:,i,ixgm));
        wp       = double(P(:,:,i,ixwm));
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
            gp       = double(P(:,:,i,ixgm));
            wp       = double(P(:,:,i,ixwm));
            cp       = double(P(:,:,i,ixcs));
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
    bp           = ((bp>th).*(slices{ixgm}+slices{ixwm}))>th;
    slices{ixgm} = slices{ixgm}.*bp;
    slices{ixwm} = slices{ixwm}.*bp;

    if niter2>0
        cp           = double(c(:,:,i));
        cp           = ((cp>th).*(slices{ixgm}+slices{ixwm}+slices{ixcs}))>th;
        slices{ixcs} = slices{ixcs}.*cp;
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4)
        P(:,:,i,k1) = slices{k1}./tot;
    end 
end