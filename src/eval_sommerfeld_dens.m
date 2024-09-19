function [u, gradu] = eval_sommerfeld_dens(zks, somm_disc, dens, targinfo)
%
%  This routine returns the Fourier transform of the densities
%  corresponding to solving the transmission problem with a collection
%  of sources and dipoles.
%
%  Input arguments:
%  -  zks: complex(2)
%      wave numbers in the upper and lower half plane
%  -  somm_disc: struct
%      struct containing details about sommerfeld discretization
%  -  dens: structure
%       structure containing sommerfeld densities
%     *  dens.mu - double layer density
%     *  dens.sig - single layer density
%  - targinfo: struct or double (2,nt)
%      structure or array containing target info
%     * targinfo.r - target locations  
%       
    if isa(targinfo, 'struct')
        targs = targinfo.r;
    elseif isa(targinfo, 'double')
        targs = targinfo;
    else
        error('EVAL_SOMMERFELD_DENS:Invalid targinfo type\n')
    end
    
    [~, nt] = size(targs);
    u = complex(zeros(nt,1));
    gradu = complex(zeros(2,nt));

    sig = dens.sig;
    mu = dens.mu;

    iind = targs(2,:) > 0;
    targ_up = targs(:,iind);
    [~, ntup] = size(targ_up);
    u_up = complex(zeros(ntup,1));
    gradu_up = complex(zeros(2,ntup));
    xfac1 = somm_disc.xfac1;
    yfac1 = somm_disc.yfac1;
    w1 = somm_disc.w1;

    
    for i=1:ntup
        zfac = exp(1j*yfac1.*targ_up(2,i)).*exp(1j*xfac1*targ_up(1,i));
        u_up(i) = sum(zfac.*w1.*(sig./(-1j*yfac1) + mu));
        gradu_up(1,i) = sum(zfac.*w1.*xfac1.*(sig./(-1j*yfac1) + mu))*1j;
        gradu_up(2,i) = sum(zfac.*w1.*yfac1.*(sig./(-1j*yfac1) + mu))*1j;
    end

    u(iind) = u_up;
    gradu(:,iind) = gradu_up;


    iind = targs(2,:) < 0;
    targs_down = targs(:,iind);
    [~, ntdown] = size(targs_down);
    u_down = complex(zeros(ntdown,1));
    gradu_down = complex(zeros(2,ntdown));
    xfac2 = somm_disc.xfac2;
    yfac2 = somm_disc.yfac2;
    w2 = somm_disc.w2;

    
    for i=1:ntdown
        zfac = exp(-1j*yfac2.*targ_down(2,i)).*exp(1j*xfac2*targ_down(1,i));
        u_down(i) = sum(zfac.*w2.*(sig./(-1j*yfac2) - mu));
        gradu_down(1,i) = sum(zfac.*w2.*xfac2.*(sig./(-1j*yfac2) - mu))*1j;
        gradu_down(2,i) = -sum(zfac.*w2.*yfac2.*(sig./(-1j*yfac2) - mu))*1j;
    end

    u(iind) = u_down;
    gradu(:,iind) = gradu_down;



end
