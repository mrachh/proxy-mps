%% Nontrivial example
n = 16;
nref = 3;

opts = [];
opts.n = n;
opts.nref = nref;
 
bl = 1.0;
xylim = [-bl/2, -bl/2; bl/2 bl/2];
zk = 10;

q = @(x,y) 50*cos(30*x).*exp(-100*(x.^2+ y.^2));

% Setup proxy points 
npxy0 = 30;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = bl*1.25/2;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr*bsize;
pw = pw*bsize;
sqpw = sqrt(pw);

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;

[Amat, L ,sols] = get_volume_scattering_matrix(zk, q, ...
      xylim, pinfo, opts);

%% Test the solution for plane wave data

dir = pi/3;
fuse = @(x,y,z) zk.^2*pw_freespace(x, y, zk, dir).*q(x,y);
usol = L.solve(fuse);

xeval = 5.5;
yeval = 4.3;
uex = usol.ext(xeval, yeval);

xx = pinfo.r(1,:).';
yy = pinfo.r(2,:).';
rnx = pinfo.n(1,:).';
rny = pinfo.n(2,:).';

uin = pw_freespace(xx, yy, zk, dir).*sqpw(:);
dudnin = 1j*zk*uin.*(cos(dir).*rnx + sin(dir).*rny);

data = [uin; dudnin];
uout = Amat*data;

upxy = uout(1:npxy).*sqpw(:);
dudnpxy = uout(npxy+1:end).*sqpw(:);

Sk = kernel('helm', 's', zk);
Dk = kernel('helm', 'd', zk);

targinfo = [];
targinfo.r = [xeval; yeval];
uval = Dk.eval(pinfo, targinfo)*upxy - Sk.eval(pinfo, targinfo)*dudnpxy;

fprintf('Error in scattering matrix soln=%d\n', abs(uval-uex));


function u = pw_freespace(x, y, zk, dir)
    u = exp(1j.*zk.*(cos(dir).*x + sin(dir).*y));
end