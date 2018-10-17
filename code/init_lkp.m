function dat = init_lkp(dat,opt)
K   = opt.template.K;
S0  = numel(dat);   
lkp = 1:K;
mg  = ones(1,K);

for s=1:S0
    dat{s}.gmm.part.lkp = lkp;
    dat{s}.gmm.part.mg  = mg;
end
%==========================================================================