%  Number of proxy points as a function of aspect ratio
%
%
clear
clc

% Estimated diameter of the object
a = 5;
b = 1/2;

nas = 3;
npxy_pts = 8;
errs = zeros(nas, npxy_pts);
npxys_all = zeros(nas, npxy_pts);

for ii=1:nas
    a = 5*2^(ii-1);
    zk = 10*pi/a;

    dsep = 1;

    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    pref = []; 
    pref.k = 16;
    chnkr = chunkerfunc(@(t) ellipse(t, a, b), cparams, pref); 
    chnkr = sort(chnkr);
    wts = chnkr.wts; wts = wts(:);

    % Set up kernel
    Ck = kernel('helm', 'c', zk, [1.0, -1j*zk]);
    Ckp = kernel('helm', 'cprime', zk, [1.0, -1j*zk]);
    
    Sk = kernel('helm', 's', zk);
    Skp = kernel('helm', 'sprime', zk);
    
    Dk = kernel('helm', 'd', zk);
    Dkp = kernel('helm', 'dprime', zk);
    
    n = chnkr.npt;
    A = chunkermat(chnkr, Ck) + 0.5*eye(n);

    % Analytic solution test
    ns = 10;
    ns_out = ns; 
    nt = 1000;
    [sources, sources_out, targets] = get_sources(a, b, ns, ns_out, nt);
    strengths = randn(ns, 1);

    err1 = test_analytic_solution(zk, chnkr, sources, targets, strengths, A);
    fprintf('error in analytic solution test=%d\n',err1);

    
    ls = [a+b, b+dsep/3];
    
    
    for jj = 1:npxy_pts
        fprintf('\n\nWavenumber = %d\n',zk)
        fprintf('npxys_id = %d\n',jj)
        % npxys = [70, 20];
        
        Nw_hor = 30*(jj+3);
        Nw_ver = 3*(jj+3);

        npxys = [Nw_hor, Nw_ver];
    

        opts = [];
        opts.iflege = 1;
        [pr,ptau,pw,pin] = chnk.flam.proxy_rect_pts(ls, npxys, opts);

        pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
    
        pinfo = [];
        pinfo.r = pr;
        pinfo.n = pn;
        pinfo.wts = pw;
        sqpw = sqrt(pw);

    
    %% Now construct the scattering matrix
    
        chnkr1 = chunkerfunc(@(t) ellipse(t, a, b), cparams, pref); 
        chnkr1 = sort(chnkr1);
    
        shift = [0, 2*b + dsep];
        chnkr2 = chunkerfunc(@(t) ellipse(t, a, b, shift), cparams, pref); 
        chnkr2 = sort(chnkr2);
    
        A11 = A;
        pinfo1 = pinfo;
        [S1, C1, Cp1] = get_scattering_matrices(zk, chnkr1, pinfo1, A11);
    
        pinfo2 = pinfo;
        pinfo2.r = pinfo2.r + shift';
    
        A22 = chunkermat(chnkr2, Ck) + 0.5*eye(n);
        [S2, C2, Cp2] = get_scattering_matrices(zk, chnkr2, pinfo2, A22);
        

        % Set up plane wave data as sources;
        thet = 0;%rand*2*pi;
        ct = cos(thet); sd = sin(thet);
    
        px1 = pinfo1.r(1,:).';
        py1 = pinfo1.r(2,:).';
    
        pnx1 = pinfo1.n(1,:).';
        pny1 = pinfo1.n(2,:).';


        pw1 = pinfo1.wts;
        pw1 = pw1(:);
        sqpw1 = sqrt(pw1);
    
        uincs1 = exp(1j * zk * (px1*ct + py1*sd) );
        gxuincs1 = 1j * zk * ct .* uincs1;
        gyuincs1 = 1j * zk * sd .* uincs1;
    
        dudnincs1 = pnx1 .* gxuincs1 + pny1 .* gyuincs1;


        px2 = pinfo2.r(1,:).';
        py2 = pinfo2.r(2,:).';
    
        pnx2 = pinfo2.n(1,:).';
        pny2 = pinfo2.n(2,:).';
    
        pw2 = pinfo2.wts;
        pw2 = pw2(:);
        sqpw2 = sqrt(pw2);
    
        uincs2 = exp(1j * zk * (px2*ct + py2*sd));
        gxuincs2 = 1j * zk * ct .* uincs2;
        gyuincs2 = 1j * zk * sd .* uincs2;
    
        dudnincs2 = pnx2 .* gxuincs2 + pny2 .* gyuincs2;

        uin_pxy = [uincs1 .* sqpw1(:); 
                  dudnincs1 .* sqpw1(:); 
                 uincs2 .* sqpw2(:); 
                dudnincs2 .* sqpw2(:)]; 
    
        [mm,~] = size(S1);
    
        Stotal = [S1, zeros(mm);
                  zeros(mm),  S2];

    % Test Stotal by solving two decoupled problems on the individual
    % scatterers
    
    uin1_pxy = [uincs1.*sqpw1(:); 
                dudnincs1.*sqpw1(:)];
       
    uin2_pxy = [uincs2.*sqpw2(:); 
                dudnincs2.*sqpw2(:)];    
    
    % Build the translation operators
    D_p1_to_p2 = sqpw.' .* Dk.eval(pinfo1, pinfo2) .* sqpw;
    S_p1_to_p2 = sqpw.' .* Sk.eval(pinfo1, pinfo2) .* sqpw; 
    Dp_p1_to_p2 = sqpw.' .* Dkp.eval(pinfo1, pinfo2) .* sqpw;
    Sp_p1_to_p2 = sqpw.' .* Skp.eval(pinfo1, pinfo2) .* sqpw;

    D_p2_to_p1 = sqpw.' .* Dk.eval(pinfo2, pinfo1) .* sqpw; 
    Dp_p2_to_p1 = sqpw.' .* Dkp.eval(pinfo2, pinfo1) .* sqpw; 
    S_p2_to_p1 = sqpw.' .* Sk.eval(pinfo2, pinfo1) .* sqpw;
    Sp_p2_to_p1 = sqpw.' .* Skp.eval(pinfo2, pinfo1) .* sqpw;
    
    
    [nnn, ~] = size(D_p1_to_p2);
    zmat = zeros(nnn);
    zeye = eye(nnn);
    Tmat = [zmat,        zmat,          D_p2_to_p1,  -S_p2_to_p1;
            zmat,        zmat,          Dp_p2_to_p1, -Sp_p2_to_p1;
            D_p1_to_p2,  -S_p1_to_p2,   zmat,        zmat;
            Dp_p1_to_p2, -Sp_p1_to_p2   zmat,        zmat];


    % Build the matrix
    
    [nn, ~] = size(Stotal);
    Ssolve = eye(nn) - (Stotal + eye(nn))*Tmat;
    
    udata_pxy = Stotal * uin_pxy;    
    uout_pxy = Ssolve \ udata_pxy;
    
    % Construct solution via direct computation
    chnkrs(1,2) = chunker();
    chnkrs(1) = chnkr1;
    chnkrs(2) = chnkr2;
    
    chnkrtotal = merge(chnkrs);
    x = chnkrtotal.r(1,:).';
    y = chnkrtotal.r(2,:).';
    ubdry = -exp(1j * zk * (x*ct + y*sd) ); 
    
    ntot = chnkrtotal.npt;
   
    Afull2 = chunkermat(chnkrtotal, Ck) + 0.5*eye(ntot);
    sig = Afull2 \ ubdry;
    
    cinfo_use = [];
    cinfo_use.r = chnkrtotal.r(:,:);
    cinfo_use.n = chnkrtotal.n(:,:);
    
    pinfo_use = [];
    npxy1 = size(px1,1);
    npxy2 = size(px2,1);
    
    pinfo_use.r = zeros(2,npxy1+npxy2);
    pinfo_use.n = zeros(2,npxy1+npxy2);
    
    pinfo_use.r(1,1:npxy1) = px1;
    pinfo_use.r(1,npxy1+1:end) = px2;
    
    pinfo_use.r(2,1:npxy1) = py1;
    pinfo_use.r(2,npxy1+1:end) = py2;
    
    pinfo_use.n(1,1:npxy1) = pnx1;
    pinfo_use.n(1,npxy1+1:end) = pnx2;
    
    pinfo_use.n(2,1:npxy1) = pny1;
    pinfo_use.n(2,npxy1+1:end) = pny2;
    
    sqpw_total = [sqpw1(:); sqpw2(:)];
    
    u_ex = sqpw_total(:).*(Ck.eval(cinfo_use, pinfo_use)*(sig.*chnkrtotal.wts(:)));
    dudn_ex = sqpw_total(:).*(Ckp.eval(cinfo_use, pinfo_use)*(sig.*chnkrtotal.wts(:)));


    uout_pxy_ex = [u_ex(1:npxy1);
                   dudn_ex(1:npxy1);
                   u_ex((npxy1+1):end);
                   dudn_ex((npxy1+1):end)];
               
    err1 = norm(uout_pxy_ex - uout_pxy);
    fprintf('Error in final solution = %d\n',err1);
    errs(ii,jj) = err1;
    npxys_all(ii,jj) = size(pinfo.r, 2);
    

    end %jj

end %ii

%% Plot the results

figure
clf
semilogy(npxys_all(1,:), errs(1,:), 'k.', 'MarkerSize', 20); hold on;
semilogy(npxys_all(1,:), errs(2,:), 'b.', 'MarkerSize', 20); 
semilogy(npxys_all(1,:), errs(3,:), 'r.', 'MarkerSize', 20); 
ylim([10^-15, 1])
xlim([200, 800])
xticks([200, 400, 600, 800])
ss = '\fontsize{12}{0}\selectfont';
xlabel('\fontsize{15}{0}\selectfont $n_{p}$', 'Interpreter','latex');
ylabel('\fontsize{15}{0}\selectfont $\varepsilon_{a}$', 'Interpreter','latex');
legend([ss '$a =10$'], [ss '$a = 20$'], [ss '$a = 40$'], 'interpreter', 'latex','Location','SouthWest')

set(gca, 'FontSize', 15)
savefig(gcf, 'a_dep');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig, 'a_dep', 'pdf')




function [S, C_chnkr_to_pxy, Cp_chnkr_to_pxy] = get_scattering_matrices(zk, chnkr, pinfo, A)
    wts = chnkr.wts(:);
    Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
    Ckp = kernel('helm', 'cprime', zk, [1.0,-1j*zk]);
    Sk = kernel('helm', 's', zk);
    Dk = kernel('helm', 'd', zk);

    pr = pinfo.r;
    pw = pinfo.wts;
    
    [~,npxy] = size(pr);
    sqpw = sqrt(pw);
    S = zeros(2*npxy);
    
    
    cinfo = [];
    cinfo.r = chnkr.r(:,:);
    cinfo.n = chnkr.n(:,:);
    
    D_pxy_to_chnkr = Dk.eval(pinfo, cinfo).*(sqpw(:).');
    S_pxy_to_chnkr = Sk.eval(pinfo, cinfo).*(sqpw(:).');


    C_chnkr_to_pxy = sqpw(:).*Ck.eval(cinfo, pinfo).*(wts(:).');
    Cp_chnkr_to_pxy = sqpw(:).*Ckp.eval(cinfo, pinfo).*(wts(:).');

    AtmpD = A \ D_pxy_to_chnkr;
    AtmpS = A \ S_pxy_to_chnkr;

    iind1 = 1:npxy;
    iind2 = (npxy+1):(2*npxy);
    S(iind1, iind1) = C_chnkr_to_pxy*AtmpD;
    S(iind1, iind2) = -C_chnkr_to_pxy*AtmpS;

    S(iind2, iind1) = Cp_chnkr_to_pxy*AtmpD;
    S(iind2, iind2) = -Cp_chnkr_to_pxy*AtmpS;
end


function [sources, sources_out, targets] = get_sources(a, b, ns, ns_out, nt)


    % sources

    
    ts = 2*pi*rand(ns, 1);
    ause = max(a,b);

    buse = min(a,b);
    sources = ellipse(ts, buse, buse);
    sources = 0.5*sources;
    

    % exterior sources
    ts = 2*pi*rand(ns_out, 1);
    sources_out = ellipse(ts, ause, ause);
    sources_out = sources_out .* (3.5 + 3*repmat(rand(1, ns_out), 2, 1));


    % targets
    ts = 2*pi*rand(nt, 1);
    targets = ellipse(ts, ause, ause);
    targets = targets .* (1 + 3*repmat(rand(1, nt), 2, 1));
end


function err1 = test_analytic_solution(zk, chnkr, sources, targets, strengths, A)
    
    Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
    Sk = kernel('helm', 's', zk);
    
   
    srcinfo  = []; srcinfo.r = sources;
    targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
    kernmats = Sk.eval(srcinfo, targinfo);
    ubdry_analytic = kernmats*strengths;

    % Solve
    sol_analytic = A \ ubdry_analytic;


    % Compute exact solution
    srcinfo  = []; srcinfo.r  = sources;
    targinfo = []; targinfo.r = targets;
    kernmatstarg = Sk.eval(srcinfo, targinfo);
    utarg = kernmatstarg*strengths;

    % Compute solution using chunkerkerneval
    % evaluate at targets and compare

    opts.usesmooth = false;
    opts.verb = false;
    opts.quadkgparams = {'RelTol', 1e-16, 'AbsTol', 1.0e-16};


    Dsol = chunkerkerneval(chnkr, Ck, sol_analytic, targets, opts);


    wchnkr = chnkr.wts;
    err1 = norm(utarg-Dsol, 'inf') / dot(abs(sol_analytic(:)), wchnkr(:));

end



function [r,d,d2] = ellipse(t,a,b, varargin)
    shift = [0 0];
    if (nargin == 4)
        shift = varargin{1};
    end
    x0 = shift(1);
    y0 = shift(2);
    xs = x0 + a*cos(t);
    ys = y0 + b*sin(t);
    dxs = -a*sin(t);
    dys = b*cos(t);

    d2xs = -a*cos(t);
    d2ys = -b*sin(t);

    r = [(xs(:)).'; (ys(:)).'];
    d = [(dxs(:)).'; (dys(:)).'];
    d2 = [(d2xs(:)).'; (d2ys(:)).'];
end