%
% Aspect ratio depedence in layered medium case
%

clear
clc


dsep = 1;




cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;

a = 5;
b = 1/2;


shift0 = [0, 2];

nas = 3;
npxy_pts = 8;
errs = zeros(nas, npxy_pts);
npxys_all = zeros(nas, npxy_pts);


for ii = 1:nas
    a = 5*2^(ii-1);
    zk = 10*pi/a;


    chnkr = chunkerfunc(@(t)ellipse(t, a, b, shift0), cparams, pref); 
    chnkr = sort(chnkr);
    wts = chnkr.wts; wts = wts(:);


    dsep = 1;
    zks = zk*[1,1.3];
    
    zk1 = zks(1);
    zk2 = zks(2);
    coefs = [-2, 2j*zk1];
    
    xylim = [-2*a, 2*a; 1.0, 4*b+5*dsep/3+2.0];
    tol = 1e-12;
    somm_disc = get_sommerfeld_disc(zks, xylim, tol);
    
    fkerns = @(s,t) eval_lm_smat(zks, somm_disc, s, t);
    fkernsp = @(s,t) eval_lm_spmat(zks, somm_disc, s, t);
    fkernd = @(s,t) eval_lm_dmat(zks, somm_disc, s, t);
    fkerndp = @(s,t) eval_lm_dpmat(zks, somm_disc, s, t);
    fkernc = @(s,t) coefs(1)*eval_lm_dmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_smat(zks, somm_disc, s, t);
    fkerncp = @(s,t) coefs(1)*eval_lm_dpmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_spmat(zks, somm_disc, s, t);

    ifquad = 1;
    tic, Amat = coefs(2)*eval_lm_mat(zks, somm_disc, chnkr, chnkr, 's', ifquad);
    Amat = Amat + coefs(1)*eval_lm_mat(zks, somm_disc, chnkr, chnkr, 'd', ifquad); toc;
    [~, na] = size(Amat);
    Amat = Amat - eye(na);

%% Now test scattering matrix for 2 objects
    
    chnkr1 = chunkerfunc(@(t) ellipse(t, a, b, shift0), cparams, pref); 
    chnkr1 = sort(chnkr1);

    Amat1 = Amat;
    
    shifty = [0, 2*b + dsep];
    shift = shift0 + shifty;
    chnkr2 = chunkerfunc(@(t) ellipse(t, a, b, shift), cparams, pref); 
    chnkr2 = sort(chnkr2);

    ifquad = 1;
    tic, Amat2 = coefs(2)*eval_lm_mat(zks, somm_disc, chnkr2, chnkr2, 's', ifquad);
    Amat2 = Amat2 + coefs(1)*eval_lm_mat(zks, somm_disc, chnkr2, chnkr2, 'd', ifquad); toc;
    [~, na] = size(Amat2);
    Amat2 = Amat2 - eye(na);



%% Construct solution via direct computation
    chnkrs(1,2) = chunker();
    chnkrs(1) = chnkr1;
    chnkrs(2) = chnkr2;
    
    chnkrtotal = merge(chnkrs);
    
    
    cinfo_use = [];
    cinfo_use.r = chnkrtotal.r(:,:);
    cinfo_use.n = chnkrtotal.n(:,:);
    cinfo_use.wts = chnkrtotal.wts(:);
    % ubdry = -eval_incident_field_halfspace(zks, somm_disc, s, cinfo_use, strengths);
    
    ifquad = 1;
    tic, Afull = coefs(2)*eval_lm_mat(zks, somm_disc, chnkrtotal, chnkrtotal, 's', ifquad);
    Afull = Afull + coefs(1)*eval_lm_mat(zks, somm_disc, chnkrtotal, chnkrtotal, 'd', ifquad); toc;
    [~, ntot] = size(Afull);
    Afull = Afull - eye(ntot);
   
    ls = [a+b, b+dsep/3];
    
    
    for jj = 1:npxy_pts
        fprintf('\n\nWavenumber = %d\n',zk)
        fprintf('npxys_id = %d\n',jj)
        
        Nw_hor = 30*(jj+3);
        Nw_ver = 3*(jj+3);

        npxys = [Nw_hor, Nw_ver];
    

        opts = [];
        opts.iflege = 1;
        [pr, ptau, pw, ~] = chnk.flam.proxy_rect_pts(ls, npxys, opts);

        pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
    
        pinfo = [];
        pinfo.r = pr + shift0.';
        pinfo.n = pn;
        pinfo.wts = pw;
        sqpw = sqrt(pw);        
        
        
        A11 = Amat1;
        pinfo1 = pinfo;
        [S1, C1, Cp1] = get_scattering_matrices_half_space(zks, somm_disc, chnkr1, pinfo1, A11);
        
        pinfo2 = pinfo;
        pinfo2.r = pinfo2.r + shifty';
        
        A22 = Amat2;
        [S2, C2, Cp2] = get_scattering_matrices_half_space(zks, somm_disc, chnkr2, pinfo2, A22);
        
        
        [mm,~] = size(S1);
        
        Stotal = [S1, zeros(mm);
                  zeros(mm),  S2];
        
        
        
        
        %% Set up boundary data
        pw1 = pinfo1.wts;
        pw1 = pw1(:);
        sqpw1 = sqrt(pw1);
        pn1 = pinfo1.n;
        
       
        thet = pi/3;
        alpha = -thet;
        [uincs1, gradu] = planewave(zks, pinfo1, alpha);
        dudn = gradu(1,:).*pn1(1,:) + gradu(2,:).*pn1(2,:);
        dudnincs1 = dudn(:);
        
        
        pw2 = pinfo2.wts;
        pw2 = pw2(:);
        sqpw2 = sqrt(pw2);
        pn2 = pinfo2.n;
        
        [uincs2, gradu] = planewave(zks, pinfo2, alpha);
        dudn = gradu(1,:).*pn2(1,:) + gradu(2,:).*pn2(2,:);
        dudnincs2 = dudn(:);
        
        
        uin_pxy = [uincs1 .* sqpw1(:); 
                  dudnincs1 .* sqpw1(:); 
                  uincs2 .* sqpw2(:); 
                  dudnincs2 .* sqpw2(:)]; 
        
        %% Build the translation operators
        
        % Build the translation operators
        D_p1_to_p2 = sqpw.' .* fkernd(pinfo1, pinfo2) .* sqpw;
        S_p1_to_p2 = sqpw.' .* fkerns(pinfo1, pinfo2) .* sqpw; 
        Dp_p1_to_p2 = sqpw.' .* fkerndp(pinfo1, pinfo2) .* sqpw;
        Sp_p1_to_p2 = sqpw.' .* fkernsp(pinfo1, pinfo2) .* sqpw;
        
        D_p2_to_p1 = sqpw.' .* fkernd(pinfo2, pinfo1) .* sqpw; 
        S_p2_to_p1 = sqpw.' .* fkerns(pinfo2, pinfo1) .* sqpw;
        Dp_p2_to_p1 = sqpw.' .* fkerndp(pinfo2, pinfo1) .* sqpw; 
        Sp_p2_to_p1 = sqpw.' .* fkernsp(pinfo2, pinfo1) .* sqpw;
          
        
        %%
        
        [nnn, ~] = size(D_p1_to_p2);
        zmat = zeros(nnn);
        zeye = eye(nnn);
        Tmat = [zmat,        zmat,          D_p2_to_p1,  -S_p2_to_p1;
                zmat,        zmat,          Dp_p2_to_p1, -Sp_p2_to_p1;
                D_p1_to_p2,  -S_p1_to_p2,   zmat,        zmat;
                Dp_p1_to_p2, -Sp_p1_to_p2   zmat,        zmat];
        
        
        %% Build the matrix
        
        [nn, ~] = size(Stotal);
        Ssolve = eye(nn) - (Stotal + eye(nn))*Tmat;
        udata_pxy = Stotal * uin_pxy;
        
        
        uout_pxy = Ssolve \ udata_pxy;
        
        ubdry = -planewave(zks, chnkrtotal, alpha);
        sig = Afull \ ubdry;
        
        
        pinfo_use = [];
        npxy = size(pinfo1.r, 2);
        npxy1 = npxy;
        npxy2 = npxy;
        
        pinfo_use.r = zeros(2,npxy1+npxy2);
        pinfo_use.n = zeros(2,npxy1+npxy2);
        
        pinfo_use.r(:,1:npxy1) = pinfo1.r(:,:);
        pinfo_use.r(:,npxy1+1:end) = pinfo2.r(:,:);
        
        
        pinfo_use.n(:,1:npxy1) = pinfo1.n;
        pinfo_use.n(:,npxy1+1:end) = pinfo2.n;
        
        
        sqpw_total = [sqpw1(:); sqpw2(:)];
        
        %%
        
        C_chnkrtot_to_pxy = sqpw_total(:).*fkernc(chnkrtotal, pinfo_use);
        Cp_chnkrtot_to_pxy = sqpw_total(:).*fkerncp(chnkrtotal, pinfo_use);
        
        sig = sig.*chnkrtotal.wts(:);
        
        u_ex = C_chnkrtot_to_pxy*sig;
        dudn_ex = Cp_chnkrtot_to_pxy*sig;
        
        
        uout_pxy_ex = [u_ex(1:npxy1);
                       dudn_ex(1:npxy1);
                       u_ex((npxy1+1):end);
                       dudn_ex((npxy1+1):end)];
        %%           
        err1 = norm(uout_pxy_ex - uout_pxy);
        fprintf('Error in final solution = %d\n',err1);
        errs(ii,jj) = err1;
        npxys_all(ii,jj) = size(pinfo.r, 2);
    end
end


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
legend([ss '$a =10$'], [ss '$a = 20$'], [ss '$a = 40$'], 'interpreter', 'latex', 'Location', 'SouthWest');

set(gca, 'FontSize', 15)
savefig(gcf, 'a_lm_dep');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig, 'a_lm_dep', 'pdf')




function [sources_out, targets] = get_sources_halfspace(narms, amp, ns, nt, shift0)


    % sources
    tmin = pi/6;
    
   
    % exterior sources
    ts = tmin + (pi-2*tmin)*rand(ns, 1);
    sources_out = starfish(ts, narms, amp);
    sources_out = sources_out .* (4 + 3*repmat(rand(1, ns), 2, 1)) + shift0.';


    % targets
    ts = tmin + (pi-2*tmin)*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (4 + 3*repmat(rand(1, nt), 2, 1)) + shift0.';
end





function [S, C_chnkr_to_pxy, Cp_chnkr_to_pxy] = get_scattering_matrices_half_space(zks, somm_disc, chnkr, pinfo, Amat)
    coefs = [-2, 2j*zks(1)];
    dk1 = kernel('helm', 'd', zks(1));
    sk1 = kernel('helm', 's', zks(1));

    pr = pinfo.r;
    pw = pinfo.wts;
    pn = pinfo.n;
    
    [~,npxy] = size(pr);
    sqpw = sqrt(pw);
    S = complex(zeros(2*npxy));
    
    cinfo = [];
    cinfo.r = chnkr.r(:,:);
    cinfo.n = chnkr.n(:,:);

    D_pxy_to_chnkr = dk1.eval(pinfo, cinfo).*(sqpw(:).');
    S_pxy_to_chnkr = sk1.eval(pinfo, cinfo).*(sqpw(:).');

    % the C and Cp matrices need to be updated to include the correction
    fkern = @(s,t) coefs(1)*eval_lm_dmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_smat(zks, somm_disc, s, t);
    fkernp = @(s,t) coefs(1)*eval_lm_dpmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_spmat(zks, somm_disc, s, t);

    C_chnkr_to_pxy = sqpw(:).*fkern(chnkr, pinfo).*chnkr.wts(:).';
    Cp_chnkr_to_pxy = sqpw(:).*fkernp(chnkr, pinfo).*chnkr.wts(:).';

    AtmpD = Amat \ D_pxy_to_chnkr;
    AtmpS = Amat \ S_pxy_to_chnkr;

    iind1 = 1:npxy;
    iind2 = (npxy+1):(2*npxy);
    S(iind1, iind1) = C_chnkr_to_pxy*AtmpD;
    S(iind1, iind2) = -C_chnkr_to_pxy*AtmpS;

    S(iind2, iind1) = Cp_chnkr_to_pxy*AtmpD;
    S(iind2, iind2) = -Cp_chnkr_to_pxy*AtmpS;


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

