function sched = get_sched(opt)
def = spm_shoot_defaults;

if opt.template.do
    sched.gmm = opt.gmm.niter;
    sched.reg = def.sched;
    sched.eul = def.eul_its(3:end);
else    
    sched.gmm = opt.gmm.niter;       
    sched.reg = def.sched;      
    sched.eul = Inf;
end

sched.a = def.sched(2:end);

sched.labels = opt.gmm.labels.S;
%==========================================================================