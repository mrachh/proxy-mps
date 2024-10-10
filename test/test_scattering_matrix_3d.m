clear
clc

zk = 4.1;
S = geometries.sphere(1.0, 3, [], 6, 11);

%%
% Set up kernel

rep_pars = [-1j*zk; 1.0];
Ck = @(s,t) helm3d.kern(zk, s, t, 'c', rep_pars);
Ckp = @(s,t) helm3d.kern(zk, s, t, 'cprime', rep_pars);
Sk = @(s,t) helm3d.kern(zk, s, t, 's');
Skp = @(s,t) helm3d.kern(zk, s, t, 'sprime');
Dk = @(s,t) helm3d.kern(zk, s, t, 'd');

eps = 1e-7;
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S, ...
   eps, zk, rep_pars, S, opts_quad);

spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
Q.spmat = spmat;
%%
zpars = [zk; rep_pars(1); rep_pars(2)];
P = zeros(S.npts,1);
A = helm3d.dirichlet.matgen(1:S.npts, 1:S.npts, S, zpars, P, Q);
%%
% Analytic solution test

% sources

ns = 10;
thet = pi*rand(ns, 1);
phi = 2*pi*rand(ns, 1);
sources = zeros(3,ns);
rs = 0.5*rand(ns, 1);
sources(1,:) = rs.*sin(thet).*cos(phi);
sources(2,:) = rs.*sin(thet).*sin(phi);
sources(3,:) = rs.*cos(thet);
strengths = randn(ns, 1);

% exterior sources
ns_out = ns;
thet = pi*rand(ns, 1);
phi = 2*pi*rand(ns, 1);
sources_out = zeros(3,ns);
rs = 1.3 + 0.5*rand(ns, 1);
sources_out(1,:) = rs.*sin(thet).*cos(phi);
sources_out(2,:) = rs.*sin(thet).*sin(phi);
sources_out(3,:) = rs.*cos(thet);


% targets
nt = 100;
thet = pi*rand(nt, 1);
phi = 2*pi*rand(nt, 1);
targets = zeros(3, nt);
rs = 1.3 + 0.5*rand(nt, 1);
targets(1,:) = rs.*sin(thet).*cos(phi);
targets(2,:) = rs.*sin(thet).*sin(phi);
targets(3,:) = rs.*cos(thet);


% Set up boundary data
srcinfo  = []; srcinfo.r = sources;
targinfo = []; targinfo.r = S.r(:,:); targinfo.n = S.n(:,:);
kernmats = Sk(srcinfo, targinfo);
ubdry_analytic = kernmats*strengths;

% Solve
sol_analytic = A \ ubdry_analytic;
sol_analytic2 = solver(S, 'helm', 'dir', ubdry_analytic, eps, zk, rep_pars);



%% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = Sk(srcinfo, targinfo);
utarg = kernmatstarg*strengths;

% Compute solution using smooth quadrature weights and compare

srcinfo = []; srcinfo.r = S.r(:,:); srcinfo.n = S.n(:,:);
targinfo = []; targinfo.r = targets;
kern_eval = Ck(srcinfo, targinfo);
sol_use = sol_analytic.*S.wts(:);
ucomp = kern_eval*(sol_use);

relerr = norm(utarg-ucomp, 'inf') / dot(abs(sol_analytic(:)), S.wts(:));
fprintf('relative l_inf/l_1 error %5.2e\n', relerr);



%% Setup proxy points 
npxy0 = 30;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = 1.25;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr*bsize;
pr = pr + shift.';
pw = pw*bsize;
sqpw = sqrt(pw);

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;




figure(1)
clf
plot(chnkr,'kx'), hold on;
plot(pr(1,:),pr(2,:),'b.');
plot(sources_out(1,:),sources_out(2,:),'ro');
axis equal

% Now start building scattering matrix


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
S(iind1, iind1) = -C_chnkr_to_pxy*AtmpD;
S(iind1, iind2) = C_chnkr_to_pxy*AtmpS;

S(iind2, iind1) = -Cp_chnkr_to_pxy*AtmpD;
S(iind2, iind2) = Cp_chnkr_to_pxy*AtmpS;


%% Now test scattering matrix for 

% Set up boundary data
srcinfo  = []; srcinfo.r = sources_out;
ubdry_pxy = (Sk.eval(srcinfo, pinfo)*strengths).*sqpw(:);
dudnbdry_pxy = (Skp.eval(srcinfo, pinfo)*strengths).*(sqpw(:));

% Verify green's theorem
ubdry = -D_pxy_to_chnkr*(ubdry_pxy) + S_pxy_to_chnkr*(dudnbdry_pxy);
uex = Sk.eval(srcinfo, cinfo)*strengths;
fprintf('error in greens id= %d\n', norm(ubdry-uex)/norm(uex));

% find projection onto plane wave data

u_data = [ubdry_pxy; dudnbdry_pxy];
u_sol = S*(u_data);

sig = A \ uex;

% Note C already has left and right weight scaling 
uex_sol = (C_chnkr_to_pxy*(sig));
dudnex_sol = (Cp_chnkr_to_pxy*(sig));

utotex = [uex_sol; dudnex_sol];

fprintf('error in scattering matrix solution=%d\n',norm(utotex-u_sol)./norm(u_sol));
