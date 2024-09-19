clear
zks = 1*pi*[1, 1.3]/2;


% Setup proxy points 
npxy0 = 50;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = 1.1;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);
pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
pr = pr * bsize;
pw = pw * bsize;

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;
sqpw = sqrt(pw);


nobjs = 2;
shifts = zeros(2,nobjs);
pinfo_list = cell(nobjs,1);

for iobj = 1:nobjs
    shifts(1,iobj) = 4*bsize*(-nobjs/2 + iobj) + rand*0.3*bsize;
    shifts(2,iobj) = 3;
    pinfo_list{iobj} = pinfo;
    pinfo_list{iobj}.r = pinfo.r + shifts(:,iobj);
end


pstruct = cat(1, pinfo_list{:});

srcuse = [];
rr = horzcat(pstruct.r);
nn = horzcat(pstruct.n);
ww = vertcat(pstruct.wts);

ptot = [];
ptot.r = rr;
ptot.n = nn;
ptot.wts = ww;

sqpwtot = sqrt(ww);


plot(rr(1,:), rr(2,:), 'k.'); hold on;

xmin = min(rr(1,:));
xmax = max(rr(1,:));

xuse = max(xmax, abs(xmin));

ymin = min(rr(2,:));
ymax = max(rr(2,:));
xylim = [-2*xuse, 2*xuse; ymin, ymax];
tol = 1e-6;
somm_disc = get_sommerfeld_disc(zks, xylim, tol);

targinfo = [];
nt = 100;
targinfo.r = rand(2,nt);
targinfo.r(1,:) = targinfo.r(1,:)*(xmax - xmin) + xmin;
targinfo.r(2,:) = targinfo.r(2,:)*(ymax - ymin) + ymin;

thet = 2*pi*rand(1,nt);
targinfo.n = zeros(2,nt);
targinfo.n(1,:) = cos(thet);
targinfo.n(2,:) = sin(thet);

plot(targinfo.r(1,:), targinfo.r(2,:), 'b.');

ntest = 100;
ttest = [];
ttest.r = targinfo.r(:,1:ntest);
ttest.n = targinfo.n(:,1:ntest);

ttest = targinfo;


%%
fkerns = @(s,t) eval_lm_smat(zks, somm_disc, s, t);
fkernsp = @(s,t) eval_lm_spmat(zks, somm_disc, s, t);
fkernd = @(s,t) eval_lm_dmat(zks, somm_disc, s, t);
fkerndp = @(s,t) eval_lm_dpmat(zks, somm_disc, s, t);

Tmat = complex(zeros(2*nt, 2*nobjs*npxy));
iind = 1:2*npxy;
for i = 1:nobjs
    istart = (i-1)*2*npxy;

    Dkmat = fkernd(pinfo_list{i}, ttest).*sqpw.';
    Skmat = fkerns(pinfo_list{i}, ttest).*sqpw.';
    Dkpmat = fkerndp(pinfo_list{i}, ttest).*sqpw';
    Skpmat = fkernsp(pinfo_list{i}, ttest).*sqpw'; 
    Tmat(:, istart+iind) = [Dkmat, -Skmat; ...
                            Dkpmat, -Skpmat];

end
Tmat(isnan(Tmat)) = 0;
%%

n = npxy;
ntot = 2*n*nobjs;
density = rand(ntot,1).*repmat(sqpw, [2*nobjs,1]);
out_ex = Tmat*density;

eps = 1e-10;
tic, [u, gradu] = eval_lm_targ_fmm(eps, zks, somm_disc, nobjs, n, pinfo_list, sqpw, targinfo, density); toc



out_test = complex(zeros(size(out_ex)));
out_test(1:ntest) = u(1:ntest);
out_test(ntest+1:end) = gradu(1,1:ntest).*ttest.n(1,:) + ...
    gradu(2,1:ntest).*ttest.n(2,:);

err = norm(out_test - out_ex)/norm(density);
fprintf('error in matvec=%d\n',err);

xmin = min(rr(1,:)); xmax = max(rr(1,:));
ymin = min(rr(2,:)); ymax = max(rr(2,:));
xlam = (xmax - xmin)*zks(2)/2/pi;
ylam = (ymax - ymin)*zks(2)/2/pi;

fprintf('Number of wavelengths in x and y direction = %d %d\n', xlam, ylam);