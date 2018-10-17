function dat = init_armijo(dat)
S0 = numel(dat);
for s=1:S0
    [~,~,~,C] = obs_info(dat{s});
    
    dat{s}.armijo.bf   = ones(1,C);
    dat{s}.armijo.aff  = 1;
    dat{s}.armijo.nl   = 1;
    dat{s}.armijo.prop = 1;
end
%==========================================================================