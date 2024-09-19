function [dens] = solve_sommerfeld_dens(zks, somm_disc, srcinfo)
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
%  -  srcinfo: structure
%        structure containing sourceinfo
%     *  srcinfo.sources: double(2,n)    
%           source locations, $x_{j}$
%     *  srcinfo.charges: complex(n,1) 
%           charge densities, $c_{j}$ (optional, 
%           default - term corresponding to charges dropped)
%     *  srcinfo.dipstr: complex(n,1)
%           dipole densities, $d_{j}$ (optional, 
%           default - term corresponding to dipoles dropped)
%     *  srcinfo.dipvec: double(2,n) 
%           dipole orientation vectors, $v_{j}$ (optional
%           default - term corresponding to dipoles dropped) 
%       
%
    sources = srcinfo.sources;
    [~, ns] = size(sources);
    if isfield(srcinfo, 'charges')
        charges = srcinfo.charges;
        charges = charges(:);
    else
        charges = complex(zeros(ns,1));
    end

    if isfield(srcinfo, 'dipvec') && isfield(srcinfo, 'dipstr')
        dipstr = srcinfo.dipstr;
        dipstr = dipstr(:);
        dipvec = srcinfo.dipvec;
    else
        dipvec = zeros(2,ns);
        dipstr = complex(zeros(ns,1));
    end
%
%  Create incident field due to sources on top
%
    iind = sources(2,:) > 0;
    sources_up = sources(:,iind);
    charges_up = charges(iind);
    dipstr_up = dipstr(iind);
    dipvec_up = dipvec(:,iind);

    xfac1 = somm_disc.xfac1;
    yfac1 = somm_disc.yfac1;

    [~, nup] = size(sources_up);
    n = somm_disc.n;
    sommer_exp1 = complex(zeros(n,1));
    for i=1:nup
        zdot = dipstr_up(i)*1j*(dipvec_up(1,i)*xfac1 - dipvec_up(2,i)*yfac1);
        zinc = exp(-1j*xfac1*sources_up(1,i)).*exp(1j*yfac1*sources_up(2,i));
        sommer_exp1 = sommer_exp1 + charges_up(i)*zinc;
        sommer_exp1 = sommer_exp1 - zinc.*zdot;
    end

    iind = sources(2,:) < 0;
    
    sources_down = sources(:,iind);
    charges_down = charges(iind);
    dipstr_down = dipstr(iind);
    dipvec_down = dipvec(:,iind);
    [~, ndown] = size(sources_down);

    xfac2 = somm_disc.xfac2;
    yfac2 = somm_disc.yfac2;
    sommer_exp2 = complex(zeros(n,1));

    for i=1:ndown
        zdot = dipstr_down(i)*1j*(dipvec_down(1,i)*xfac2 + dipvec_down(2,i)*yfac2);
        zinc = exp(-1j*xfac2*sources_down(1,i)).*exp(-1j*yfac2*sources_down(2,i));
        sommer_exp2 = sommer_exp2 + charges_down(i)*zinc;
        sommer_exp2 = sommer_exp2 - zinc.*zdot;
    end

    sq1 = -1j*yfac1;
    sq2 = -1j*yfac2;
    ff = -sommer_exp1./sq1 + sommer_exp2./sq2;
    gg = -(sommer_exp1 + sommer_exp2);
    det = -2 - (sq2./sq1 + sq1./sq2);
    mu = (-2*ff + (1./sq2 - 1./sq1).*gg)./det;
    sig = ((sq1 - sq2).*ff + 2*gg)./det;
    
    dens = [];
    dens.mu = mu;
    dens.sig = sig;

end
