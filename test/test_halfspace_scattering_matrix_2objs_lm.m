clear
clc


% Estimated diameter of the object
diam = 2;

% desired size of object in wavelengths
dwav = 2; 

% separation as measured in wavelenths
dsep = 1;

zk = dwav*2*pi/diam;

zks = zk*[1,1.3];

zk1 = zks(1);
zk2 = zks(2);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;

shift0 = [0.3, 2.4];
chnkr = chunkerfunc(@(t) starfish(t, narms, amp, shift0), cparams, pref); 
chnkr = sort(chnkr);
wts = chnkr.wts; wts = wts(:);

coefs = [-2, 2j*zk1];
xylim = [-15, 15; 0.5, 10];
tol = 1e-9;
somm_disc = get_sommerfeld_disc(zks, xylim, tol);

fkerns = @(s,t) eval_lm_smat(zks, somm_disc, s, t);
fkernsp = @(s,t) eval_lm_spmat(zks, somm_disc, s, t);
fkernd = @(s,t) eval_lm_dmat(zks, somm_disc, s, t);
fkerndp = @(s,t) eval_lm_dpmat(zks, somm_disc, s, t);
fkernc = @(s,t) coefs(1)*eval_lm_dmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_smat(zks, somm_disc, s, t);
fkerncp = @(s,t) coefs(1)*eval_lm_dpmat(zks, somm_disc, s, t) + coefs(2)*eval_lm_spmat(zks, somm_disc, s, t);

Amat = chunkermat(chnkr, fkernc);
[~, na] = size(Amat);
Amat = Amat - eye(na);

% Analytic solution test
ns = 10;
nt = 10;
[sources_out, targets] = get_sources_halfspace(narms, amp, ns, nt, shift0);
strengths = randn(ns, 1);


% Setup proxy points 
npxy0 = 50;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = 1.1;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr * bsize + shift0.';
pw = pw * bsize;

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;
sqpw = sqrt(pw);

% plot boundary, proxy and sources
figure(1)
clf
plot(chnkr,'kx'), hold on;
plot(pr(1,:),pr(2,:),'b.');
plot(sources_out(1,:),sources_out(2,:),'ro');
axis equal


%% Now test scattering matrix for 2 objects

chnkr1 = chunkerfunc(@(t) starfish(t, narms, amp, shift0), cparams, pref); 
chnkr1 = sort(chnkr1);

shiftx = [5 0];
shift = shiftx + shift0;
chnkr2 = chunkerfunc(@(t) starfish(t, narms, amp, shift), cparams, pref); 
chnkr2 = sort(chnkr2);


A11 = Amat;
pinfo1 = pinfo;
[S1, C1, Cp1] = get_scattering_matrices_half_space(zks, somm_disc, chnkr1, pinfo1, A11);

pinfo2 = pinfo;
pinfo2.r = pinfo2.r + shiftx';

A22 = Amat;
[S2, C2, Cp2] = get_scattering_matrices_half_space(zks, somm_disc, chnkr2, pinfo2, A22);


[mm,~] = size(S1);

Stotal = [S1, zeros(mm);
          zeros(mm),  S2];


figure(2)
clf
plot(chnkr1,'k.'); hold on;
plot(chnkr2, 'r.');
plot(pinfo1.r(1,:),pinfo1.r(2,:),'ko');
plot(pinfo2.r(1,:),pinfo2.r(2,:),'ro');
plot(sources_out(1,:), sources_out(2,:),'b.');

axis equal

%% Set up boundary data
pw1 = pinfo1.wts;
pw1 = pw1(:);
sqpw1 = sqrt(pw1);
pn1 = pinfo1.n;

s = []; s.r = sources_out;

[uincs1, gradu] = eval_incident_field_halfspace(zks, somm_disc, s, pinfo1, strengths);
dudn = gradu(1,:).*pn1(1,:) + gradu(2,:).*pn1(2,:);
dudnincs1 = dudn(:);


pw2 = pinfo2.wts;
pw2 = pw2(:);
sqpw2 = sqrt(pw2);
pn2 = pinfo2.n;

[uincs2, gradu] = eval_incident_field_halfspace(zks, somm_disc, s, pinfo2, strengths);
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


%% Test Green's identity
s2 = [];
s2.r = [5.0; 3.1];
[vin, gradv] = eval_incident_field_halfspace(zks, somm_disc, s2, pinfo2);
dvindn = gradv(1,:).*pn2(1,:) + gradv(2,:).*pn2(2,:);
dvindn = dvindn(:);

vex = eval_incident_field_halfspace(zks, somm_disc, s2, pinfo1);
vex = vex.*sqpw1(:);

vc = D_p2_to_p1*(vin.*sqpw2(:)) - S_p2_to_p1*(dvindn.*sqpw2(:));
fprintf('error in translation operator=%d\n', norm(vc - vex));


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
afun = @(x)  Ssolve*x;
uout_pxy_gmres = gmres(afun, udata_pxy, [], 1e-14);

%% Construct solution via direct computation
chnkrs(1,2) = chunker();
chnkrs(1) = chnkr1;
chnkrs(2) = chnkr2;

chnkrtotal = merge(chnkrs);


cinfo_use = [];
cinfo_use.r = chnkrtotal.r(:,:);
cinfo_use.n = chnkrtotal.n(:,:);
cinfo_use.wts = chnkrtotal.wts(:);
ubdry = -eval_incident_field_halfspace(zks, somm_disc, s, cinfo_use, strengths);

Afull = chunkermat(chnkrtotal, fkernc);
ntot = chnkrtotal.npt;
Afull = Afull - eye(ntot);

sig = Afull \ ubdry;


pinfo_use = [];
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

