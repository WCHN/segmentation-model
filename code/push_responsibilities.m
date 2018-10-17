function [Z_a,bb] = push_responsibilities(Z,y,dm_a,mat_a,bb)
K = size(Z,4);
for k=1:K
    [Z_k,c] = spm_diffeo('push',Z(:,:,:,k),y,dm_a(1:3));
    
    if k==1
        % Calculate bounding box from count image (c)
        [bb,xbb,ybb,zbb,dm_bb] = get_bb(c,mat_a,dm_a,bb);                
        Z_a                    = zeros([dm_bb K],'single');
    end       

    Z_a(:,:,:,k) = Z_k(xbb,ybb,zbb);
end
%==========================================================================