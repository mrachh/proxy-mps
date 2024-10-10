function [Amat, varargout] = get_volume_scattering_matrix(zk, qfun, ...
      xylim, pinfo, opts)
%
%  This function constructs the scattering matrix for a spatially
%  varying sound speed. The sound speed is assumed to compactly
%  supported within xylim where we use a heirarchical Poincare Steklov
%  solver
%
%  Input arguments:
%  -  zks: complex
%      wave number in free space
%  -  qfun: function handle 
%      qfun(x,y) is the function handle for the spatially varying 
%      sound speed
%  -  xylim: double (2,2)
%      [xmin, xmax; ymin, ymax]
%      xy extent of the spatially varying sound speed
%  - pinfo: struct
%      pinfo.r - location of proxy points
%      pinfo.n - normals on the proxy surface
%      pinfo.wts - weights for integrating smooth functions
%        on the proxy surface
%  - opts: struct
%      structure with optional input arguments
%      opts.n - order of discretization (16)
%      opts.nref - number of refinements of base patch,
%           default value determined based on points per wavelength
%           considerations, min value 2, assumes
%           qfun is resolved by such a mesh.
%
%  Output arguments:
%   - Amat: complex(2*np, 2*np)
%       scattering matrix
%   - L: surfacebie object
%       compressed operator for evaluating solutions for new
%       data
%
%
    if nargin < 5
        opts = [];
    end

    n = 16;
    if isfield(opts, 'n')
        n = opts.n;
    end

    xlam = (xylim(2,1) - xylim(1,1))*real(zk)/2/pi;
    ylam = (xylim(2,2) - xylim(1,2))*real(zk)/2/pi;
    lam = max(xlam, ylam);
    nref = max(ceil(log(lam)/log(2)) + 1, 2);
    if isfield(opts, 'nref')
        nref = opts.nref;
    end
    
    quse = @(x,y,z) qfun(x, y);

    x = @(u,v) u;
    y = @(u,v) v;
    z = @(u,v) 0*u;
    rect = xylim(:).';
    

    L = surfacebie(n, x, y, z, q=quse, zk=zk, nref=nref, rect=rect);

    varargout{1} = L;

    [~ ,np] = size(pinfo.r);
    Amat = complex(zeros(2*np,2*np));
    xeval = pinfo.r(1,:).';
    yeval = pinfo.r(2,:).';
    Dkp = kernel('helm', 'dprime', zk);
    Skp = kernel('helm', 'sprime', zk);

    opts_quad = [];
    opts_quad.accel = true;
    
    
    Dkp_u  = @(f) chunkerkerneval(L.chnkr, Dkp, f, pinfo, opts_quad);
    Skp_du = @(f) chunkerkerneval(L.chnkr, Skp, f, pinfo, opts_quad);
    sqpw = sqrt(pinfo.wts(:));

    sols = cell(2*np,1);
    for i = 1:np
        fuse = @(x,y,z) -kern_eval(x, y, zk, pinfo.r(:,i), ...
            pinfo.n(:,i), pinfo.wts(i), 'd', qfun);
        sols{i} = L.solve(fuse);
        uval = sols{i}.ext(xeval, yeval).*sqpw;
        dudnval = Dkp_u(sols{i}.dir) - Skp_du(sols{i}.neu);
        dudnval = dudnval.*sqpw;
        Amat(1:np,i) = uval;
        Amat(np+1:end,i) = dudnval;
    end

    for i=1:np
        fuse = @(x,y,z) kern_eval(x, y, zk, pinfo.r(:,i), ...
            pinfo.n(:,i), pinfo.wts(i), 's', qfun);
        sols{i+np} = L.solve(fuse);
        uval = sols{i+np}.ext(xeval, yeval).*sqpw;
        dudnval = Dkp_u(sols{i+np}.dir) - Skp_du(sols{i+np}.neu);
        dudnval = dudnval.*sqpw;
        Amat(1:np,i+np) = uval;
        Amat(np+1:end,i+np) = dudnval;
    end
    varargout{2} = sols;
end


function f = kern_eval(x, y, zk, r, rn, w, type, qfun)
    srcinfo = [];
    srcinfo.r = r;
    srcinfo.n = rn;
    xx = x(:).'; yy = y(:).';

    targinfo = [];
    targinfo.r = [xx; yy];
    
    feval = kernel('helm', type, zk);
    f = feval.eval(srcinfo, targinfo);
    wuse = reshape(w, [size(r,2),1]);
    f = f.*sqrt(wuse);
    f = squeeze(reshape(f, [size(x), size(r,2)]));
    qvals = qfun(x,y);
    f=  zk.^2*f.*qvals;
end