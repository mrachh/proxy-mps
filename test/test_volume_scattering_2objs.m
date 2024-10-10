clear
clc
%% Setup problem
n = 16;
nref = 3;

opts_scat = [];
opts_scat.n = n;
opts_scat.nref = nref;
 
bl = 1.0;
xylim = [-bl/2, -bl/2; bl/2 bl/2];
zk = 3;

qfun1 = @(x,y,c) 50*cos(30*(x-c(1))).*exp(-100*((x-c(1)).^2+ (y-c(2)).^2));

% Setup proxy points 
npxy0 = 30;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = bl*1.25/2;
[pr, ptau, pw, ~] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr*bsize;
pw = pw*bsize;
sqpw = sqrt(pw);

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;

%% Construct scattering matrix for 1 object

qfunuse = @(x,y) qfun1(x,y, [0,0]);
[Amat, L ,sols] = get_volume_scattering_matrix(zk, qfunuse, ...
      xylim, pinfo, opts_scat);

%%


% Set up kernel
Sk = kernel('helm', 's', zk);
Skp = kernel('helm', 'sprime', zk);

Dk = kernel('helm', 'd', zk);
Dkp = kernel('helm', 'dprime', zk);

%% Now test scattering matrix for 2 objects

pinfo1 = pinfo;
shift = [3*bl, 1.5*bl];
pinfo2 = pinfo;
pinfo2.r = pinfo2.r + shift';

xylim2 = [-bl/2 + shift(1), -bl/2 + shift(2); bl/2 + shift(1), bl/2+shift(2)];

[Amat2, L2 ,sols2] = get_volume_scattering_matrix(zk, @(x,y) qfun1(x,y,shift), ...
      xylim2, pinfo2, opts_scat);


% Set up plane wave data as sources;
thet = pi/3;
cd = cos(thet); sd = sin(thet);

px1 = pinfo1.r(1,:).';
py1 = pinfo1.r(2,:).';

pnx1 = pinfo1.n(1,:).';
pny1 = pinfo1.n(2,:).';
pw1 = pinfo1.wts;
pw1 = pw1(:);
sqpw1 = sqrt(pw1);

uincs1 = exp(1j*zk*(px1*cd + py1*sd));
gxuincs1 = 1j*zk*cd.*uincs1;
gyuincs1 = 1j*zk*sd.*uincs1;

dudnincs1 = pnx1.*gxuincs1 + pny1.*gyuincs1;


px2 = pinfo2.r(1,:).';
py2 = pinfo2.r(2,:).';

pnx2 = pinfo2.n(1,:).';
pny2 = pinfo2.n(2,:).';
pw2 = pinfo2.wts;
pw2 = pw2(:);
sqpw2 = sqrt(pw2);

uincs2 = exp(1j*zk*(px2*cd + py2*sd));
gxuincs2 = 1j * zk * cd .* uincs2;
gyuincs2 = 1j * zk * sd .* uincs2;

dudnincs2 = pnx2 .* gxuincs2 + pny2 .* gyuincs2;

uin_pxy = [uincs1.*sqpw1; 
          dudnincs1.*sqpw1(:); 
          uincs2.*sqpw2(:); 
          dudnincs2.*sqpw2(:)]; 

[mm,~] = size(Amat);

Stotal = [Amat, zeros(mm);
          zeros(mm),  Amat];


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
%%
[nn, ~] = size(Stotal);

udata_pxy = Stotal*uin_pxy;

afun = @(x) x - (Stotal + eye(nn))*(Tmat*x);
uout_pxy_gmres = gmres(afun, udata_pxy, [], 1e-14);

%% Construct solution via direct computation

qfun_tot = @(x,y,z) qfun1(x,y,[0,0]) + qfun1(x,y,shift);
rect = [-bl/2, -bl/2; bl/2+shift(1), bl/2 + shift(2)];

nref = 4;
x = @(u,v) u;
y = @(u,v) v;
z = @(u,v) 0*u;
rect = rect(:).';
    
Ltot = surfacebie(n, x, y, z, q=qfun_tot, zk=zk, nref=nref, rect=rect);

fuse = @(x,y,z) zk.^2.*qfun_tot(x,y,z).*exp(1j.*zk.*(x.*cd + y.*sd));
utot = Ltot.solve(fuse);

%%
rext = [-3.1*bl/2; -5.1*bl/2];
uex = utot.ext(rext(1), rext(2));

upxy1 = uout_pxy_gmres(1:npxy).*sqpw1(:);
upxy2 = uout_pxy_gmres(2*npxy+1:3*npxy).*sqpw2(:);

dudnpxy1 = uout_pxy_gmres(npxy+1:2*npxy).*sqpw1(:);
dudnpxy2 = uout_pxy_gmres(3*npxy+1:end).*sqpw2(:);

targinfo = [];
targinfo.r = rext;
u1 = Dk.eval(pinfo1, targinfo)*upxy1 - Sk.eval(pinfo1, targinfo)*dudnpxy1;
u2 = Dk.eval(pinfo2, targinfo)*upxy2 - Sk.eval(pinfo2, targinfo)*dudnpxy2;

u = u1 + u2;
fprintf('error in 2 obj scattering=%d\n', abs(u-uex));

%% Plot things
figure(1)
clf
plot(pinfo1.r(1,:), pinfo1.r(2,:), 'k.'); hold on
plot(pinfo2.r(1,:), pinfo2.r(2,:), 'b.'); hold off


q = surfacefun(qfun_tot, utot.int.domain);
figure(2)
clf
plot(q)



