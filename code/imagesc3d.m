function imagesc3d(img,varargin)
% FORMAT imagesc3d(img,varargint)
% img      - A 2D or 3D image
% varargin - Options for MATLAB's imagesc() function
%
% Same as imagesc() if 2D, otherwise displays a montage of slices from a 3D
% image.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

dm = size(img);
dm = [dm 1]; 

if numel(dm) > 3 && dm(4) > 1
    error('imagesc3d does not support >3D arrays!')
elseif dm(3) > 1
    
    % Montage parameters: this is the spacing for picking slices in the
    % z-axis
    Spacing = 8;
    
    % Set up montage
    z  = 1:Spacing:dm(3);
    N  = numel(z);    
    nr = floor(sqrt(N));
    nc = ceil(N/nr);  
    
    % Create montage
    mtg = zeros([nr*dm(1) nc*dm(2)],'single');
    cnt = 1;
    for r=1:nr
        for c=1:nc
            if cnt > numel(z)
                break
            end
            
            mtg(1 + (r - 1)*dm(1):r*dm(1),1 + (c - 1)*dm(2):c*dm(2)) = img(:,:,z(cnt));
            
            cnt = cnt + 1;
        end
    end   
    img = mtg;
    clear mtg
end

% Show image with imagesc()
imagesc(img,varargin{:});
%==========================================================================