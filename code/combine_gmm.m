function Z = combine_gmm(Z,part)
dm  = size(Z);
I   = dm(1);
K   = dm(2);
K_p = numel(part);
if K==K_p
    return
end

nZ = zeros([I K],'single');
for k=1:K
    nZ(:,k) = sum(Z(:,part==k),2);
end  
Z = nZ;
%==========================================================================