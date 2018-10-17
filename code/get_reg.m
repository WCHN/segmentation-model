function reg = get_reg(dat,fct,opt)

% FFT of Green's function
reg.prm                             = [fct.subj.vs fct.subj.ff*opt.reg.rparam];
if opt.reg.int_args > 1, reg.Greens = spm_shoot_greens('kernel',fct.subj.dm(1:3),reg.prm);
else,                    reg.Greens = [];
end

% Get initial velocities
if isnumeric(dat.reg.v)
    % Initial velocity stored in array
    reg.v = dat.reg.v;                            
else
    % Initial velocity read from NIfTI
    reg.v = single(dat.reg.v.dat(:,:,:,:));       
end

% Make deformation
reg.y      = make_deformation(reg.v,reg.prm,opt.reg.int_args,reg.Greens);

% Affine matrix
E          = spm_dexpm(dat.reg.r,opt.reg.B);
reg.Affine = fct.templ.mat\E*fct.subj.mat;

%==========================================================================