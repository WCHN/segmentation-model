function dat = meancorrect_aff(dat,opt)
% FORMAT dat = meancorrect_aff(dat,opt)
% dat - Subjects data structure
% opt - Options structure
%
% Zero-centre affine parameters across subjects.
% (Using: dat.reg.r, opt.reg.mc_maff, opt.reg.B)

if opt.reg.mc_aff            
    S0 = numel(dat);
    
    r = zeros(size(opt.reg.B,3),S0);
    for s=1:S0
        r(:,s) = dat{s}.reg.r;
    end

    r_avg = mean(r,2);
    clear r

    for s=1:S0
        dat{s}.reg.r = dat{s}.reg.r - r_avg;
    end            
end
%==========================================================================