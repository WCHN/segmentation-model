function [subsmp,o,dm] = get_subsampling_grid(dm,vs,samp)

if samp == 0
    subsmp    = [];
%     subsmp.MT = eye(4);
    o         = 1;
    return
end

sk  = max([1 1 1],round(samp*[1 1 1]./vs)); 

% For multiplying and dividing displacements to map from the subsampled voxel indices
% and the actual image voxel indices.
sk4 = reshape(sk,[1 1 1 3]);

% Mapping from indices of subsampled voxels to indices of voxels in image(s).
MT = [sk(1) 0 0 (1-sk(1));0 sk(2) 0 (1-sk(2)); 0 0 sk(3) (1-sk(3));0 0 0 1];

subsmp.sk  = sk;
subsmp.sk4 = sk4;
subsmp.MT  = MT;

[x0,y0,o] = ndgrid(1:sk(1):dm(1),1:sk(2):dm(2),1);
z0        = 1:sk(3):dm(3);

subsmp.x0 = x0;
subsmp.y0 = y0;
subsmp.z0 = z0;    

dm        = [size(subsmp.x0) length(subsmp.z0)]; % Sub-sampled dimensions
subsmp.dm = dm;