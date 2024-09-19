clear
clc

zk = 10.1;

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
shift = [1500 0];
start = tic; chnkr = chunkerfunc(@(t) starfish(t, narms, amp, shift), cparams, pref); 
t1 = toc(start);
chnkr = sort(chnkr);
wts = chnkr.wts; wts = wts(:);

% Set up kernel
Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
Ckp = kernel('helm', 'cprime', zk, [1.0,-1j*zk]);
Sk = kernel('helm', 's', zk);
Skp = kernel('helm', 'sprime', zk);
Dk = kernel('helm', 'd', zk);

n = chnkr.npt;
A = chunkermat(chnkr, Ck) + 0.5*eye(n);

% Analytic solution test

% sources

ns = 10;
ts = 2*pi*rand(ns, 1);
sources = starfish(ts, narms, amp);
sources = 0.5*sources + shift.';
strengths = randn(ns, 1);

% exterior sources
ns_out = ns;
ts = 2*pi*rand(ns_out, 1);
sources_out = starfish(ts, narms, amp);
sources_out = sources_out .* (3 + 3*repmat(rand(1, ns_out), 2, 1)) + shift.';


% targets
nt = 1000;
ts = 2*pi*rand(nt, 1);
targets = starfish(ts, narms, amp);
targets = targets .* (1 + 3*repmat(rand(1, nt), 2, 1)) + shift.';


% Set up boundary data
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


start = tic;
Dsol = chunkerkerneval(chnkr, Ck, sol_analytic, targets, opts);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)


wchnkr = chnkr.wts;
relerr  = norm(utarg-Dsol) / (sqrt(chnkr.nch)*norm(utarg));
relerr2 = norm(utarg-Dsol, 'inf') / dot(abs(sol_analytic(:)), wchnkr(:));
fprintf('relative frobenius error %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error %5.2e\n', relerr2);



% Setup proxy points 
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
