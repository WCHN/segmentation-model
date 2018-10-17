function [bb,xbb,ybb,zbb,dm_bb] = get_bb(img,mat,dm,bb)
if nargin<4, bb = zeros(3,2); end

if ~isempty(dm) && dm(3) == 1
    bb    = [1 dm(1); 1 dm(2); 1 1];    
    xbb   = bb(1,1):bb(1,2);
    ybb   = bb(2,1):bb(2,2);
    zbb   = bb(3,1):bb(3,2);
    dm_bb = [numel(xbb) numel(ybb) numel(zbb)];

    return
end

if sum(bb(:))==0
    bb = spm_imbasics('compute_bb',img,mat,dm(1:3),'nz');

    mn = [bb(1,:) 1]';
    mx = [bb(2,:) 1]';

    mn = mat\mn;
    mx = mat\mx;

    bb = [mn(1:3), mx(1:3)];
    bb = sort(bb,2);
end

xbb   = bb(1,1):bb(1,2);
ybb   = bb(2,1):bb(2,2);
zbb   = bb(3,1):bb(3,2);
dm_bb = [numel(xbb) numel(ybb) numel(zbb)];
%==========================================================================