function [opt,holly] = default_opt(opt)

if nargin < 1, opt = struct; end

START_NL_TEMPL = 1;
def            = spm_shoot_defaults;

% holly
if ~isfield(opt,'holly') 
    opt.holly = struct;
end
holly = distribute_default(opt.holly);

% opt
if ~isfield(opt,'dir_output')
    opt.dir_output = './output/';
end

% opt.model
if ~isfield(opt,'model') 
    opt.model         = struct;
end
if ~isfield(opt.model,'tol') 
    opt.model.tol     = 1e-4;
end
if ~isfield(opt.model,'niter') 
    opt.model.niter   = 1;
end
if ~isfield(opt.model,'nam_cls') 
    opt.model.nam_cls = {};
end
if ~isfield(opt.model,'clean_up') 
    opt.model.clean_up = true;
end

% opt.model.PropPrior
if ~isfield(opt.model,'PropPrior') 
    opt.model.PropPrior    = struct;
end
if ~isfield(opt.model.PropPrior,'do') 
    opt.model.PropPrior.do = true;
end

% opt.gmm
if ~isfield(opt,'gmm') 
    opt.gmm                = struct;
end
if ~isfield(opt.gmm,'niter') 
    opt.gmm.niter          = 20;
end
if ~isfield(opt.gmm,'tol') 
    opt.gmm.tol            = 1e-4;
end
if ~isfield(opt.gmm,'pth_GaussPrior') 
    opt.gmm.pth_GaussPrior = '';
end
if ~isfield(opt.gmm,'pth_PropPrior') 
    opt.gmm.pth_PropPrior  = '';
end

% opt.gmm.hist
if ~isfield(opt.gmm,'hist') 
    opt.gmm.hist             = struct;
end
if ~isfield(opt.gmm.hist,'niter_main') 
    opt.gmm.hist.niter_main  = 0; % 6
end
if ~isfield(opt.gmm.hist,'niter_gmm') 
    opt.gmm.hist.niter_gmm   = 10; 
end
if ~isfield(opt.gmm.hist,'init_ix') 
    % For setting indices of classes (e.g. to match different modalities)
    map                      = containers.Map;
    opt.gmm.hist.init_ix     = map; 
end
if ~isfield(opt.gmm.hist,'verbose') 
    opt.gmm.hist.verbose     = true; % [true,false]
end
if ~isfield(opt.gmm.hist,'verbose_gmm') 
    opt.gmm.hist.verbose_gmm = 0; % [0,1,2,3]
end

% opt.gmm.GaussPrior
if ~isfield(opt.gmm,'GaussPrior') 
    opt.gmm.GaussPrior             = struct;
end
if ~isfield(opt.gmm.GaussPrior,'constrained') 
    opt.gmm.GaussPrior.constrained = true;
end
if ~isfield(opt.gmm.GaussPrior,'verbose') 
    opt.gmm.GaussPrior.verbose     = true; % [true,false]
end

% opt.gmm.labels
if ~isfield(opt.gmm,'labels') 
    opt.gmm.labels     = struct;
end
if ~isfield(opt.gmm.labels,'cm') 
    % For including labels (if provided)
    map                = containers.Map;
    opt.gmm.labels.cm  = map;
end
if ~isfield(opt.gmm.labels,'use') 
    opt.gmm.labels.use = false;
end
if ~isfield(opt.gmm.labels,'S') 
    opt.gmm.labels.S   = 0.98;
end

% opt.sched
if ~isfield(opt,'sched') 
    opt.sched = get_sched(opt,START_NL_TEMPL);
end

% opt.reg
if ~isfield(opt,'reg') 
    opt.reg           = struct;
end
if ~isfield(opt.reg,'rparam0') 
    opt.reg.rparam0   = 1e-1*[0 0.001 0.5 0.05 0.2];
end
if ~isfield(opt.reg,'rparam') 
    opt.reg.rparam    = opt.reg.rparam0;
    opt.reg.rparam(3) = opt.sched.reg(1)*opt.reg.rparam(3);
end
if ~isfield(opt.reg,'int_args') 
    opt.reg.int_args  = opt.sched.eul(1);
end
if ~isfield(opt.reg,'niter') 
    opt.reg.niter     = 3;
end
if ~isfield(opt.reg,'tol') 
    opt.reg.tol       = 1e-4;
end
if ~isfield(opt.reg,'strt_nl') 
    opt.reg.strt_nl   = START_NL_TEMPL;
end
if ~isfield(opt.reg,'mc_aff') 
    opt.reg.mc_aff    = false;
end
if ~isfield(opt.reg,'aff_type') 
    opt.reg.aff_type  = 'affine'; % ['translation','rotation','rigid','similitude','affine']
end
if ~isfield(opt.reg,'aff_reg') 
    opt.reg.aff_reg   = 1e0;
end
if ~isfield(opt.reg,'do_aff') 
    opt.reg.do_aff    = true;
end
if ~isfield(opt.reg,'do_nl') 
    opt.reg.do_nl     = true;
end

% opt.template
if ~isfield(opt,'template') 
    opt.template              = struct;
end
if ~isfield(opt.template,'do')
    opt.template.do           = false;
end
if ~isfield(opt.template,'pth_template')
    opt.template.pth_template = '';
end
if ~isfield(opt.template,'K')
    if exist(opt.template.pth_template,'file')
        Nii                   = nifti(opt.template.pth_template);
        opt.template.K        = Nii.dat.dim(4);        
    else
        opt.template.K        = 6;
    end
end
if ~isfield(opt.template,'vs')
    opt.template.vs           = 1.5;
end
if numel(opt.template.vs) == 1
    opt.template.vs           = opt.template.vs*ones(1,3);
end
if ~isfield(opt.template,'reg0')
    opt.template.reg0         = def.sparam;
end
if ~isfield(opt.template,'reg')
    opt.template.reg          = opt.sched.a(1)*opt.template.reg0;
end
if ~isfield(opt.template,'crop')
    opt.template.crop         = false;
end
if ~isfield(opt.template,'load_a_der')
    opt.template.load_a_der   = true;
end
if ~isfield(opt.template,'R')
    opt.template.R            = null(ones(1,opt.template.K));
end
if ~isfield(opt.template,'rem_neck')
    opt.template.rem_neck     = true;
end
if ~isfield(opt.template,'sym')
    opt.template.sym          = 'all';
end
if ~isfield(opt.template,'les_cl')
    opt.template.les_cl       = 6;
end
if ~isfield(opt.template,'niter')
    opt.template.niter        = 16;
end
if ~isfield(opt.template,'verbose')
    opt.template.verbose      = 0; % [0,1,2]
end

% opt.seg
if ~isfield(opt,'seg') 
    opt.seg                = struct;
end
if ~isfield(opt.seg,'niter')
    opt.seg.niter          = 20;
end
if ~isfield(opt.seg,'tol')
    opt.seg.tol            = 1e-4;
end
if ~isfield(opt.seg,'show')
    opt.seg.show           = false;
end
if ~isfield(opt.seg,'samp')
    opt.seg.samp           = 0;
end
if ~isfield(opt.seg,'write_mllabels')
    opt.seg.write_mllabels = true;
end

% opt.seg.mrf
if ~isfield(opt.seg,'mrf') 
    opt.seg.mrf          = struct;
end
if ~isfield(opt.seg.mrf,'ml')
    opt.seg.mrf.ml       = true;
end
if ~isfield(opt.seg.mrf,'val_diag')
    opt.seg.mrf.val_diag = 0.5;
end
if ~isfield(opt.seg.mrf,'alpha')
    opt.seg.mrf.alpha    = 1e5;
end

% opt.bf
if ~isfield(opt,'bf') 
    opt.bf               = struct;
end
if ~isfield(opt.bf,'biasfwhm')
    opt.bf.biasfwhm      = 60;
end
if ~isfield(opt.bf,'niter')
    opt.bf.niter         = 8;
end
if ~isfield(opt.bf,'tol')
    opt.bf.tol           = 1e-4;
end
if ~isfield(opt.bf,'mc_bf')
    opt.bf.mc_bf         = false;
end
if ~isfield(opt.bf,'biasreg')
    opt.bf.biasreg       = 1e5;
end
if ~isfield(opt.bf,'do')
    opt.bf.do            = true;
end
if ~isfield(opt.bf,'mc_bf_verbose')
    opt.bf.mc_bf_verbose = false; % [true,false]
end

% opt.prop
if ~isfield(opt,'prop') 
    opt.prop       = struct;
end
if ~isfield(opt.prop,'niter')
    opt.prop.niter = 1;
end
if ~isfield(opt.prop,'tol')
    opt.prop.tol   = 1e-4;
end
if ~isfield(opt.prop,'reg')
    opt.prop.reg   = 0;
end
if ~isfield(opt.prop,'do')
    opt.prop.do    = true;
end

% opt.nline_search
if ~isfield(opt,'nline_search') 
    opt.nline_search      = struct;
end
if ~isfield(opt.nline_search,'bf')
    opt.nline_search.bf   = 4;
end
if ~isfield(opt.nline_search,'aff')
    opt.nline_search.aff  = 4;
end
if ~isfield(opt.nline_search,'nl')
    opt.nline_search.nl   = 4;
end
if ~isfield(opt.nline_search,'prop')
    opt.nline_search.prop = 4;
end

% opt.clean
if ~isfield(opt,'clean') 
    opt.clean          = struct;
end
if ~isfield(opt.clean,'eyeballs')
    opt.clean.eyeballs = false;
end

% opt.clean.mrf
if ~isfield(opt.clean,'mrf') 
    opt.clean.mrf          = struct;
end
if ~isfield(opt.clean.mrf,'do')
    opt.clean.mrf.do       = false;
end
if ~isfield(opt.clean.mrf,'strength')
    opt.clean.mrf.strength = 2;
end
if ~isfield(opt.clean.mrf,'niter')
    opt.clean.mrf.niter    = 10;
end

% opt.start_it
if ~isfield(opt,'start_it') 
    opt.start_it            = struct;
end
if ~isfield(opt.start_it,'do_mg')
    opt.start_it.do_mg      = 1;
end
if ~isfield(opt.start_it,'do_prop')
    opt.start_it.do_prop    = 1;
end
if ~isfield(opt.start_it,'do_upd_mrf')
    opt.start_it.do_upd_mrf = 1;
end

% opt.do
if ~isfield(opt,'do') 
    opt.do            = struct;
end
if ~isfield(opt.do,'mg')
    opt.do.mg         = true;
end
if ~isfield(opt.do,'update_mrf')
    opt.do.update_mrf = true;
end
if ~isfield(opt.do,'mrf')
    opt.do.mrf        = false;
end

% opt.ct
if ~isfield(opt,'ct') 
    opt.ct            = struct;
end
if ~isfield(opt.ct,'GaussPrior')
    opt.ct.GaussPrior = [];
end
if isempty(opt.ct.GaussPrior)
    lb_prW         = struct;                                        
    lb_prW.KL_qVpV = 0;
    lb_prW.ElnDetV = zeros(1,opt.template.K);

    b  = 1e4*ones(1,opt.template.K);
    n  = 1e0*ones(1,opt.template.K);    
    MU = reshape([linspace(0,50,opt.template.K - 3) 800 -100 -1000],[1 opt.template.K]);
    W  = ones([1 1 opt.template.K]);

    opt.ct.GaussPrior = {MU,b,W,n,'CT',lb_prW,1:opt.template.K};   
end

% opt.dict
if ~isfield(opt,'dict') 
    opt.dict          = struct;
end
if ~isfield(opt.dict,'lkp')
    % For using multiple Gaussians per tissue
    map                = containers.Map;
    opt.dict.lkp       = map;
end
if ~isfield(opt.dict,'prop_excl')
    % For constraining a class to have zero resps by setting its proportion to
    % a large negative value
    map                = containers.Map;
    opt.dict.prop_excl = map;
end

% opt.verbose
if ~isfield(opt,'verbose') 
    opt.verbose       = struct;
end
if ~isfield(opt.verbose,'level') 
    opt.verbose.level = 2;
end
if ~isfield(opt.verbose,'model') 
    opt.verbose.model = 0; % [0,1,2,3]
end
if opt.verbose.level == 2
    opt.verbose.gmm   = 4; % [0,1,2,3,4,5]
    opt.verbose.reg   = 3; % [0,1,2,3]
    opt.verbose.bf    = 3; % [0,1,2,3]   
    opt.verbose.prop  = 3; % [0,1,2,3]   
    opt.verbose.mrf   = 3; % [0,1,2,3]  
elseif opt.verbose.level == 1
    opt.verbose.gmm   = 1;
    opt.verbose.reg   = 1;
    opt.verbose.bf    = 1;    
    opt.verbose.prop  = 1;
    opt.verbose.mrf   = 1;
else
    opt.verbose.gmm   = 0;
    opt.verbose.reg   = 0;
    opt.verbose.bf    = 0;    
    opt.verbose.prop  = 0;
    opt.verbose.mrf   = 0;    
end
if ~strcmpi(holly.mode,'for')
    opt.verbose.gmm  = 0;
    opt.verbose.reg  = 0;
    opt.verbose.bf   = 0;
    opt.verbose.prop = 0;
    opt.verbose.mrf  = 0;
end

% options when training
if opt.template.do 
    opt.seg.niter   = 1;
    opt.model.niter = 30;                        
    
    opt.reg.rparam0   = def.rparam;
    opt.reg.rparam    = opt.reg.rparam0;
    opt.reg.rparam(3) = opt.sched.reg(1)*opt.reg.rparam(3);

    opt.reg.strt_nl         = START_NL_TEMPL;       
    opt.start_it.do_prop    = START_NL_TEMPL;    
    opt.start_it.do_mg      = START_NL_TEMPL + 1;
    opt.start_it.do_upd_mrf = START_NL_TEMPL + 1;
    
    opt.gmm.labels.use = true;           
    opt.bf.mc_bf       = true;    
                    
    if opt.template.load_a_der
        opt.template.verbose = 1;        
    else        
        opt.template.verbose = 2;        
    end
    opt.verbose.model    = 3; 
    opt.bf.mc_bf_verbose = true;
       
    opt.dir_output = fullfile(opt.dir_output,'train');
else
    opt.dir_output = fullfile(opt.dir_output,'segment');
end
%==========================================================================

%==========================================================================
function sched = get_sched(opt,START_NL_TEMPL)
def = spm_shoot_defaults;

if opt.template.do
    sched.gmm = opt.gmm.niter;
    sched.reg = def.sched(2:end);    
    sched.eul = def.eul_its;
else    
    sched.gmm = opt.gmm.niter;    
    sched.reg = (2.^fliplr(0:9));        
    sched.eul = 9;
end

% mx      = 30;
% stps    = 20;
% sched.a = mx*exp(-(linspace(0,-log(1/mx),stps)));
sched.a = def.sched([2*ones(1,START_NL_TEMPL) 2:numel(def.sched)]);
%==========================================================================