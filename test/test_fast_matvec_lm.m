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


nobjs = 10;
shifts = zeros(2,nobjs);
pinfo_list = cell(nobjs,1);

nobtest = 2;
Tmat = complex(zeros(2*nobtest*npxy,2*nobjs*npxy));
for iobj = 1:nobjs
    shifts(1,iobj) = 4*bsize*(-nobjs/2 + iobj) + rand*0.3*bsize;
    shifts(2,iobj) = 3;
    pinfo_list{iobj} = pinfo;
    pinfo_list{iobj}.r = pinfo.r + shifts(:,iobj);
end



pstruct = cat(1, pinfo_list{:});

srcuse = [];
rr = horzcat(pstruct.r);
plot(rr(1,:), rr(2,:), 'k.');

xmin = min(rr(1,:));
xmax = max(rr(1,:));
xuse = max(xmax, abs(xmin));
ymin = min(rr(2,:));
ymax = max(rr(2,:));
xylim = [-2*xuse, 2*xuse; ymin, ymax];
tol = 1e-6;
somm_disc = get_sommerfeld_disc(zks, xylim, tol);


%%
fkerns = @(s,t) eval_lm_smat(zks, somm_disc, s, t);
fkernsp = @(s,t) eval_lm_spmat(zks, somm_disc, s, t);
fkernd = @(s,t) eval_lm_dmat(zks, somm_disc, s, t);
fkerndp = @(s,t) eval_lm_dpmat(zks, somm_disc, s, t);

iind = 1:2*npxy;
for i = 1:nobtest
    istart = (i-1)*2*npxy;
    for j=1:nobjs
        jstart = (j-1)*2*npxy;
        if (i ~= j)
            Dkmat = sqpw' .* fkernd(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Skmat = sqpw' .* fkerns(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Dkpmat = sqpw' .* fkerndp(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Skpmat = sqpw' .* fkernsp(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Tmat(istart + iind, jstart + iind) = [Dkmat, -Skmat; ...
                                                  Dkpmat, -Skpmat];
        end
    end
end

Tmat(isnan(Tmat)) = 0;
%%

Dkmat = sqpw'.*fkernd(pinfo_list{1}, pinfo_list{1}).*sqpw;
Skmat = sqpw'.*fkerns(pinfo_list{1}, pinfo_list{1}).*sqpw;
Dkpmat = sqpw'.*fkerndp(pinfo_list{1}, pinfo_list{1}).*sqpw;
Skpmat = sqpw'.*fkernsp(pinfo_list{1}, pinfo_list{1}).*sqpw;


Tself = [Dkmat, -Skmat; ...
         Dkpmat, -Skpmat];
Tself(isnan(Tself)) = 0;
%%

n = npxy;
ntot = 2*n*nobjs;
density = rand(ntot,1).*repmat(sqpw, [2*nobjs,1]);

eps = 1e-10;
tic, out = Tmat_fast_matvec_lm(eps, zks, somm_disc, nobjs, n, pinfo_list, sqpw, Tself, density); toc

out_ex = Tmat*density;

ntest = length(out_ex);
err = norm(out(1:ntest) - out_ex)/norm(density);
fprintf('error in matvec=%d\n',err);

xmin = min(rr(1,:)); xmax = max(rr(1,:));
ymin = min(rr(2,:)); ymax = max(rr(2,:));
xlam = (xmax - xmin)*zks(2)/2/pi;
ylam = (ymax - ymin)*zks(2)/2/pi;

fprintf('Number of wavelengths in x and y direction = %d %d\n', xlam, ylam);