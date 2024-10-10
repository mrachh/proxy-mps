classdef surfacebie

    properties
        dom
        op
        chnkr

        zk
        eta
        q

        DtN
        ItN
        dA
    end

    properties ( Hidden )
        P
        Q
        pre
        Skern
        Dkern
    end

    methods

        function L = surfacebie(n, xfun, yfun, zfun, opts)

            arguments
                n
                xfun = @(u,v) u
                yfun = @(u,v) v
                zfun = @(u,v) 0*u
                opts.nref = 0
                opts.rect = [-1 1 -1 1]
                opts.zk = 0
                opts.eta = []
                opts.q = @(x,y,z) 0*x
            end

            if ( isempty(opts.eta) )
                if ( opts.zk == 0 )
                    opts.eta = 1;
                else
                    opts.eta = real(opts.zk);
                end
            end

            L.zk  = opts.zk;
            L.eta = opts.eta;
            L.q   = opts.q;
            
            % Build a uniformly refined flat square and remap it according to the
            % given functions
            uv = surfacemesh.square(n, opts.nref, opts.rect);
            u = uv.x;
            v = uv.y;
            x = cell(size(u));
            y = cell(size(u));
            z = cell(size(u));
            for k = 1:length(uv)
                x{k} = xfun(u{k},v{k});
                y{k} = yfun(u{k},v{k});
                z{k} = zfun(u{k},v{k});
            end
            L.dom = surfacemesh(x, y, z);

            pdo = [];
            pdo.lap = 1;
            pdo.c = @(x,y,z) L.zk^2*(1-L.q(x,y,z));
            L.op = surfaceop(L.dom, pdo, method='ItI', eta=L.eta);
            L.op.build();
            xyz = L.op.patches{1}.xyz;

            [L.P, L.Q] = interpolators(xyz, n, opts.rect, opts.nref);

            nskel = n-2;
            L.chnkr = squarechunker(nskel, opts.nref, opts.rect);

            if ( L.zk == 0 )
                L.Skern = kernel('laplace', 's');
                L.Dkern = kernel('laplace', 'd');
            else
                L.Skern = kernel('helmholtz', 's', L.zk);
                L.Dkern = kernel('helmholtz', 'd', L.zk);
            end
            S = chunkermat(L.chnkr, L.Skern);
            D = chunkermat(L.chnkr, L.Dkern);
            I = eye(L.chnkr.npt);

            L.DtN = L.P * L.op.DtN() * L.Q;
            L.ItN = L.P / (I - L.op.ItI(L.eta));
            L.dA = decomposition(I/2 - D + S*L.DtN);
            L.pre = S;
        end

        function u = solve(L, f)

            u = [];
            % assign and evaluate the particular solution with 0
            % boundary data
            L.op.rhs = f;
            
            imp_part = L.op.patches{1}.du_part;
            neu_part = L.ItN * imp_part;

            rhs = L.pre * -neu_part;
            dir = L.dA \ rhs;
            neu = L.DtN * dir + neu_part;

            % Construct impedance boundary data
            bc = L.Q * (neu + 1i*L.eta*dir);
            u.int = L.op.solve(bc);

            opts = [];
            opts.accel = true;
            D_u  = @(x,y) reshape(chunkerkerneval(L.chnkr, L.Dkern, dir, [x(:) y(:)].', opts), size(x));
            S_du = @(x,y) reshape(chunkerkerneval(L.chnkr, L.Skern, neu, [x(:) y(:)].', opts), size(x));
            u.ext = @(x,y) D_u(x,y) - S_du(x,y);
            u.dir = dir;
            u.neu = neu;

        end

    end

end

function [P, Q] = interpolators(xyz, n, rect, nref)

    % Collect indices corresponding to each side of the outer boundary
    leftIdx  = find(abs(xyz(:,1) - rect(1)) < 1e-10);
    rightIdx = find(abs(xyz(:,1) - rect(2)) < 1e-10);
    downIdx  = find(abs(xyz(:,2) - rect(3)) < 1e-10);
    upIdx    = find(abs(xyz(:,2) - rect(4)) < 1e-10);

    % Reorder the HPS skeleton to match chunkie's panel ordering
    nskel = n-2;
    nside = nskel * 2^nref;
    P = speye(4*nside);
    P(:, [downIdx; rightIdx; flip(upIdx); flip(leftIdx)]) = P;

    % Interpolate from Chebyshev panels (HPS) to Gauss panels (chunkie)
    [xc, ~, vc] = chebpts(nskel, 1);
    [xl, ~, vl] = legpts(nskel);
    C2L = barymat(xl, xc, vc); C2L = repmat({C2L}, 4*2^nref, 1); C2L = matlab.internal.math.blkdiag(C2L{:});
    L2C = barymat(xc, xl, vl); L2C = repmat({L2C}, 4*2^nref, 1); L2C = matlab.internal.math.blkdiag(L2C{:});
    Q = L2C * P.';
    P = P * C2L;

end

function chnkr = squarechunker(n, nref, rect)

    if ( nargin < 2 )
        nref = 0;
    end

    if ( nargin < 3 )
        rect = [-1 1 -1 1];
    end

    m = 2^nref+1;
    [xa, xb, ya, yb] = dealm(rect);
    xab = linspace(xa, xb, m).';
    yab = linspace(ya, yb, m).';

    south = [ xab repmat(ya, m, 1) ];
    east  = [ repmat(xb, m, 1) yab ];
    north = [ flip(xab) repmat(yb, m, 1) ];
    west  = [ repmat(xa, m, 1) flip(yab) ];

    verts = [ south(1:m-1,:) ; east(1:m-1,:) ; north(1:m-1,:) ; west(1:m-1,:) ];
    verts = [ verts ; verts(1,:) ];

    nch = size(verts,1) - 1;

    pref = [];
    pref.k = n;
    chnkr = chunker(pref);
    chnkr = chnkr.addchunk(nch);

    t = legpts(n, [0 1]);

    for k = 1:nch
        chnkr.r(1,:,k) = (1-t)*verts(k,1) + t*verts(k+1,1);
        chnkr.r(2,:,k) = (1-t)*verts(k,2) + t*verts(k+1,2);
        len = sqrt((verts(k,1)-verts(k+1,1)).^2 + (verts(k,2)-verts(k+1,2)).^2);
        vx = (verts(k+1,1)-verts(k,1));
        vy = (verts(k+1,2)-verts(k,2));
        chnkr.d(1,:,k) = vx/2*ones(1, n);
        chnkr.d(2,:,k) = vy/2*ones(1, n);
        chnkr.d2(1,:,k) = zeros(1, n);
        chnkr.d2(2,:,k) = zeros(1, n);
        chnkr.adj(1,k) = k-1;
        chnkr.adj(2,k) = k+1;
    end

    chnkr.adj(1,1) = nch;
    chnkr.adj(2,nch) = 1;

    chnkr.n = normals(chnkr);
    chnkr.wts = weights(chnkr);

    chnkr.nstor(:,:,1:nch) = normals(chnkr);
    chnkr.wtsstor(:,1:nch) = weights(chnkr);

end
