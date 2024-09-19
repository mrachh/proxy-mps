clear
zk = 5*pi;


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

Sk = kernel('helm', 's', zk);
Dk = kernel('helm', 'd', zk);
Skp = kernel('helm', 'sp', zk);
Dkp = kernel('helm', 'dp', zk);

nlat = 5;
nobjs = nlat*nlat;
shifts = zeros(2,nobjs);
pinfo_list = cell(nobjs,1);

nobtest = 2;
Tmat = complex(zeros(2*nobtest*npxy,2*nobjs*npxy));
for i = 1:nlat
    for j = 1:nlat
        iobj = (j-1)*nlat + i;
        shifts(1,iobj) = 4*bsize*i + rand*0.3*bsize;
        shifts(2,iobj) = 4*bsize*j + rand*0.5*bsize;
        pinfo_list{iobj} = pinfo;
        pinfo_list{iobj}.r = pinfo.r + shifts(:,iobj);
    end
end

pstruct = cat(1, pinfo_list{:});

srcuse = [];
rr = horzcat(pstruct.r);
plot(rr(1,:), rr(2,:), 'k.');


iind = 1:2*npxy;
for i = 1:nobtest
    istart = (i-1)*2*npxy;
    for j=1:nobjs
        jstart = (j-1)*2*npxy;
        if (i ~= j)
            Dkmat = sqpw' .* Dk.eval(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Skmat = sqpw' .* Sk.eval(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Dkpmat = sqpw' .* Dkp.eval(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Skpmat = sqpw' .* Skp.eval(pinfo_list{j}, pinfo_list{i}) .* sqpw;
            Tmat(istart + iind, jstart + iind) = [Dkmat, -Skmat; ...
                                                  Dkpmat, -Skpmat];
        end
    end
end

Tmat(isnan(Tmat)) = 0;


Dkmat = sqpw'.*Dk.eval(pinfo, pinfo).*sqpw;
Skmat = sqpw'.*Sk.eval(pinfo, pinfo).*sqpw;
Dkpmat = sqpw'.*Dkp.eval(pinfo, pinfo).*sqpw;
Skpmat = sqpw'.*Skp.eval(pinfo, pinfo).*sqpw;


Tself = [Dkmat, -Skmat; ...
         Dkpmat, -Skpmat];
Tself(isnan(Tself)) = 0;
%%

n = npxy;
ntot = 2*n*nobjs;
density = rand(ntot,1).*repmat(sqpw, [2*nobjs,1]);
eps = 1e-7;
tic, out = Tmat_fast_matvec(eps, zk, nobjs, n, pinfo_list, sqpw, Tself, density); toc

out_ex = Tmat*density;

ntest = length(out_ex);
err = norm(out(1:ntest) - out_ex)/norm(out_ex);
fprintf('error in matvec=%d\n',err);

xmin = min(rr(1,:)); xmax = max(rr(1,:));
ymin = min(rr(2,:)); ymax = max(rr(2,:));
xlam = (xmax - xmin)*zk/2/pi;
ylam = (ymax - ymin)*zk/2/pi;

fprintf('Number of wavelengths in x and y direction = %d %d\n', xlam, ylam);