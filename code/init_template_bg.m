function model = init_template_bg(model)
dm                   = model.template.nii.dat.dim;
K                    = dm(4);
model.template.bg    = cell(1,2);
model.template.bg{1} = zeros(K,1);
model.template.bg{2} = zeros(K,1);
for k=1:K
    if dm(3) == 1
        model.template.bg{1}(k) = log(1/K);
        model.template.bg{2}(k) = log(1/K);    
    else
        model.template.bg{1}(k) = mean(mean(model.template.nii.dat(:,:,1,k)));
        model.template.bg{2}(k) = mean(mean(model.template.nii.dat(:,:,end,k)));    
    end
end
% figure;imagesc(model.template.nii.dat(:,:,1,7))
%==========================================================================