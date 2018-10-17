function show_segmentations(dat)

figname = '(SPM) Segmentations';

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);  
clf(f);

populations = spm_json_manager('get_populations',dat);
P           = numel(populations);
S0          = numel(dat);

nrows = min(S0,10); 

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
        
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)

            nii = nifti(dat{s}.pth.seg2d);    
            Z   = single(nii.dat(:,:,:,:));
            
            img = [];
            for k=1:size(Z,4)
                img = [img Z(:,:,:,k)'];
            end

            subplot(nrows,1,cnt_plots);
            imagesc(img,[0 1]); axis off xy image;
               
            cnt_plots = cnt_plots + 1;
                        
            if cnt == floor(10/P) || cnt == S0
                break
            end
            
            cnt = cnt + 1;
        end
    end
end

drawnow;
%==========================================================================