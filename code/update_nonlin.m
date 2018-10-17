function [dat,y,v,template,Greens,gain,opt] = update_nonlin(dat,model,obs,template,bf,labels,v,y,Affine,Greens,it_seg,it_template,scl,miss,opt)

if it_seg > 1
    % Non-linear registration parameters changes with iteration number (single-subject)
    opt = modify_opt(opt,it_seg);     
end

if it_template > 1
    % Non-linear registration parameters changes with iteration number (over populations)
    opt = modify_opt(opt,it_template);     
end

% Parse input
prop     = dat.gmm.prop;
cluster  = dat.gmm.cluster;
armijo   = dat.armijo.nl;
int_args = opt.reg.int_args;
verbose  = opt.verbose.reg;
part    = dat.gmm.part;

% Parameters
[dm,~,vs] = obs_info(dat);
ff        = get_ff(vs);      
prm       = [vs ff*opt.reg.rparam];
K         = size(template,2);
lkp       = [1 4 5; 4 2 6; 5 6 3];
const     = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
ix_tiny   = get_par('ix_tiny',dat.population,part.lkp,opt);

if armijo < 1e-6
    % Already found optimal solution
    dat.armijo.nl = min(armijo*1.25,1);
    gain          = 0;
    
    return; 
end

if int_args > 1
    % Large deformation
    Greens   = spm_shoot_greens('kernel',dm(1:3),prm); 
%     y        = make_deformation(v,prm,int_args,Greens);
%     Template = warp_template(model,y,Affine);    
end     

%--------------------------------------------------------------------------
% Compute objective function and its first and second derivatives
%--------------------------------------------------------------------------
    
% Init
g  = zeros([dm(1:3),3],'single');
H  = zeros([dm(1:3),6],'single');
y1 = spm_warps('transform',Affine,y);
if int_args > 1
    J  = spm_diffeo('def2jac',y1);
else
    Jz = reshape(Affine(1:3,1:3),[1 1 3 3]);
end

[dlb,mom] = gmm_img('init_lb_and_mom',miss);

% Neighborhood part
lnPzN = gmm_mrf('apply',dat.mrf);

Template0 = single(model.template.nii.dat(:,:,:,:));
for z=1:dm(3)

    % Get slices
    [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);

    if dat.mrf.do && numel(lnPzN) > K           
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end
    
    % Compute responsibilities and lb
    [Z,dlb,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);

    % Compute sufficient statistics 
    mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);

    % Map cluster responsibilities to tissue responsibilities
    Z = cluster2template(Z,part);              

    % Get gradient and Hessian of multinomial objective function
    Z = reshape(Z,[dm(1:2) 1 K]);
    
    if dat.mrf.do        
        dat.mrf.oZ(:,:,z,:) = uint8((2^8)*Z);
    end
    
    [gz,Hz] = mnom_objfun_slice(Z,Template0,y1,z,prop);

    % Compute affine gradient and Hessian
    if int_args > 1
        Jz = reshape(J(:,:,z,:,:),[dm(1:2) 3 3]);
    end

    % Rotate gradients, such that g1 = J'*g0;
    for d1=1:3
        tmp = 0;
        for d2=1:3
            tmp = tmp + Jz(:,:,d2,d1).*gz(:,:,d2);
        end
        g(:,:,z,d1) = tmp;
    end

    % Rotate Hessian, such that H2 = J'*H0*J
    % First do H1 = J'*H0
    RH  = zeros([dm(1:2),3,3],'single');
    for d1=1:3
        for d3=1:3
            tmp = 0;
            for d2=1:3
                tmp = tmp + Jz(:,:,d2,d1).*Hz(:,:,lkp(d2,d3));
            end
            RH(:,:,d1,d3) = tmp;
        end
    end

    % Then do H2 = H1*J
    for d1=1:3
        for d3=d1:3 % Only need compute an upper or lower triangle
            tmp = 0;
            for d2=1:3
                tmp = tmp + RH(:,:,d1,d2).*Jz(:,:,d2,d3);
            end
            H(:,:,z,lkp(d1,d3)) = tmp;
        end
    end
end
clear J y1 Z Template0

lbX = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, cluster{1}, cluster{2}, miss.L, mom.SS2b);    

dat.lb.X(end + 1)       = lbX;     
dat.lb.Z(end + 1)       = dlb.Z;   
if numel(part.mg) > numel(prop)
    dat.lb.mg(end + 1)  = dlb.mg;   
end
if ~isempty(labels)
    dat.lb.lab(end + 1) = dlb.lab;
end
if dat.mrf.do   
    dat.lb.ZN(end + 1)  = gmm_mrf('lowerbound',dat.mrf);
end

dat.lb = check_convergence('nl',dat.lb,1,verbose);

% Add regularisation and compute GN step
u      = spm_diffeo('vel2mom',v,prm);
g      = g + u;
Update = spm_diffeo('fmg',H, g,[prm 2 2]);
clear H g u

% Linesearch
%--------------------------------------------------------------------------
otemplate    = template;
ov           = v;
nline_search = opt.nline_search.nl;
for line_search=1:nline_search
    
    % Update initial velocity
    v = v - armijo*Update;
    if dm(3) == 1
        v(:,:,:,3) = 0;
    end
    
    % Compute new deformation
    [y,eul_its] = make_deformation(v,prm,int_args,Greens);
    
    % Compute new lower bound
    u     = spm_diffeo('vel2mom',v,prm);
    v_reg = 0;
    for i=1:3
        v_reg = v_reg + (-0.5*sum(sum(sum(sum(double(u(:,:,:,i)).*double(v(:,:,:,i)))))));
    end
    clear u
    nlb = v_reg;
    
    nlb = nlb + dat.lb.X(end) + dat.lb.lab(end) + dat.lb.lnDetbf(end) ...
          + dat.lb.MU(end) + dat.lb.A(end) + sum(dat.lb.bf_reg(end,:)) + dat.lb.aff_reg(end) + dat.lb.prop_reg(end) + dat.lb.mg(end) + dat.lb.ZN(end);          
    
    template = warp_template(model,y,Affine);    
    
    [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,scl,template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'Template',otemplate});                      
    nlb             = nlb + dlb.Z;

    % Check new lower bound
    if nlb > dat.lb.last
        armijo = min(armijo*1.25,1);
        
        dat.lb.Z(end + 1)     = dlb.Z;                 
        dat.lb.v_reg(end + 1) = v_reg;                                                

        break;
    else
        armijo = armijo*0.5;
        v      = ov;
        
        if line_search == nline_search            
            y        = make_deformation(v,prm,int_args,Greens);                        
            template = otemplate;            
        end
    end        
end
clear oTemplate ov

[dat.lb,gain] = check_convergence('nl',dat.lb,1,verbose,armijo,eul_its);    

if verbose >= 3
    show_bf_and_ivel(obs,dm,v);
end

% Set output
dat.armijo.nl   = armijo;
%==========================================================================