%
%  Test correction evaluator
%
%
eps = 1E-9;
zks = pi*[1,1.3];
zk1 = zks(1);
zk2 = zks(2);

close all


cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
shift = [0.3, 2.4];

chnkr = chunkerfunc(@(t) starfish(t, narms, amp, shift), cparams, pref); 


% sources

ns = 10;
tmin = pi/6;
ts = tmin + (pi-2*tmin)*rand(ns, 1);
sources = starfish(ts, narms, amp);
sources = 0.5*sources + shift.';
strengths = randn(ns, 1);

% exterior sources
ns_out = ns;
ts = tmin + (pi - 2*tmin)*rand(ns_out, 1);
sources_out = starfish(ts, narms, amp);
sources_out = sources_out .* (4 + 3*repmat(rand(1, ns_out), 2, 1)) + shift.';

nt = ns;
ts = tmin + (pi - 2*tmin)*rand(ns_out, 1);
targs = starfish(ts, narms, amp);
targs = targs .* (4 + 3*repmat(rand(1, ns_out), 2, 1)) + shift.';

targ_info = [];
targ_info.r = targs;


s = [];
s.r = sources_out;


figure(1)
clf
plot(chnkr,'k.'), hold on;
plot(sources_out(1,:),sources_out(2,:),'ro');
plot(targs(1,:),targs(2,:),'gx');
axis equal




%%



dk1 = kernel('helm', 'd', zk1);
sk1 = kernel('helm', 's', zk1);
skp1 = kernel('helm', 'sp', zk1);
dkp1 = kernel('helm', 'dp', zk1);

coefs = [-2, 2j*zk1];

Ck1 = kernel('helm', 'c', zk1, coefs);
Ck1p = kernel('helm', 'cp', zk1, coefs);

Amat = chunkermat(chnkr, Ck1);
n = chnkr.npt;
Amat = Amat - eye(n);

src_info = [];
src_info.r = chnkr.r(:,:);
src_info.charges = coefs(2)*chnkr.wts(:);
src_info.dipvec = chnkr.n(:,:);
src_info.dipstr = coefs(1)*chnkr.wts(:);

xylim = [-10, 10; 1.0, 10];
tol = 1e-9;
somm_disc = get_sommerfeld_disc(zks, xylim, tol);
Amat_corr = eval_sommerfeld_correction(zks, somm_disc, src_info, src_info);

Amat = Amat + Amat_corr;

rhs_use = -eval_incident_field_halfspace(zks, somm_disc, s, src_info);
mu_lm = Amat \ rhs_use;


opts = [];
opts.forcesmooth = true;
opts.accel = false;


p1 = chunkerkerneval(chnkr, Ck1, mu_lm, targ_info, opts);
p2 = eval_sommerfeld_correction(zks, somm_disc, src_info, targ_info)*mu_lm;

uscex = p1 + p2;

%% Now build and test scattering matrix


% Setup proxy points 
npxy0 = 30;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = 1.2;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr*bsize;
pr = pr + shift.';
pw = pw*bsize;
sqpw = sqrt(pw);

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;

plot(pr(1,:),pr(2,:),'b.');


%% Now start building scattering matrix


S = zeros(2*npxy);

cinfo = [];
cinfo.r = chnkr.r(:,:);
cinfo.n = chnkr.n(:,:);
D_pxy_to_chnkr = dk1.eval(pinfo, cinfo).*(sqpw(:).');
S_pxy_to_chnkr = sk1.eval(pinfo, cinfo).*(sqpw(:).');

% Construct Amat2 differently
fkern = @(s,t) coefs(1)*eval_lm_dmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_smat(zks, somm_disc, s, t);
fkernp = @(s,t) coefs(1)*eval_lm_dpmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_spmat(zks, somm_disc, s, t);
tic, Amat2 = chunkermat(chnkr, fkern); toc;
Amat2 = Amat2 - eye(n);


% the C and Cp matrices need to be updated to include the correction

src_info = [];
src_info.r = chnkr.r(:,:);
src_info.charges = coefs(2)*chnkr.wts(:);
src_info.dipvec = chnkr.n(:,:);
src_info.dipstr = coefs(1)*chnkr.wts(:);

[~, ns] = size(src_info.r);
[umat, gradumat] = eval_sommerfeld_correction(zks, somm_disc, src_info, pinfo);

nxtarg = repmat((pn(1,:)).',1,ns);
nytarg = repmat((pn(2,:)).',1,ns);
ckp_corr_mat = gradumat(:,:,1).*nxtarg + gradumat(:,:,2).*nytarg;
%%
C_chnkr_to_pxy = sqpw(:).*umat;
Cp_chnkr_to_pxy = sqpw(:).*ckp_corr_mat;

C_chnkr_to_pxy3 = sqpw(:).*fkern(chnkr, pinfo).*chnkr.wts(:).';
Cp_chnkr_to_pxy3 = sqpw(:).*fkernp(chnkr, pinfo).*chnkr.wts(:).';

AtmpD = Amat \ D_pxy_to_chnkr;
AtmpS = Amat \ S_pxy_to_chnkr;

iind1 = 1:npxy;
iind2 = (npxy+1):(2*npxy);
S(iind1, iind1) = C_chnkr_to_pxy3*AtmpD;
S(iind1, iind2) = -C_chnkr_to_pxy3*AtmpS;

S(iind2, iind1) = Cp_chnkr_to_pxy3*AtmpD;
S(iind2, iind2) = -Cp_chnkr_to_pxy3*AtmpS;


%% Now test scattering matrix for 

% Set up boundary data
srcinfo  = []; srcinfo.r = sources_out;

[u, gradu] = eval_incident_field_halfspace(zks, somm_disc, s, pinfo, strengths);
ubdry_pxy = u.*sqpw(:);
dudn = gradu(1,:).*pn(1,:) + gradu(2,:).*pn(2,:);
dudn = dudn(:);
dudnbdry_pxy = dudn.*(sqpw(:));

% Verify green's theorem
ubdry = -D_pxy_to_chnkr*(ubdry_pxy) + S_pxy_to_chnkr*(dudnbdry_pxy);
uex = eval_incident_field_halfspace(zks, somm_disc, s, cinfo, strengths);
fprintf('error in greens id= %d\n', norm(ubdry-uex)/norm(uex));

%% find projection onto plane wave data

u_data = [ubdry_pxy; dudnbdry_pxy];
u_sol = S*(u_data);

% Note that the scattering matrix corresponded to data uin instead of -uin
% earlier. To fix that, there is an extra negative sign.
sig = Amat \ (-uex);   

% Note C already has left and right weight scaling 
uex_sol = (C_chnkr_to_pxy3*(sig));
dudnex_sol = (Cp_chnkr_to_pxy3*(sig));

utotex = [uex_sol; dudnex_sol];

fprintf('error in scattering matrix solution=%d\n',norm(utotex-u_sol)./norm(u_sol));

%% Test subroutinized version of constructing scattering matrices

[S2, C_chnkr_to_pxy2, Cp_chnkr_to_pxy2] =  ...
    get_scattering_matrices_half_space(zks, somm_disc, chnkr, pinfo, Amat);

return

function [f,fd,fdd] = flat_interface(t, a, b, t0, t1)
    
    phi   = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid  = @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd = @(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f = zeros([2,size(t)]);
    fd = zeros([2,size(t)]);
    fdd = zeros([2,size(t)]);
    
    f(1,:) = t + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 
    fd(1,:)= 1 + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) = 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));
    
    f(2,:) = 0;
    fd(2,:) = 0;
    fdd(2,:) = 0;
        
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

