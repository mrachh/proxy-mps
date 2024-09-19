function out = Tmat_fast_matvec_lm(eps, zks, somm_disc, npanels, nnn, plist, sqpw, Tself, density)
% Tmat_fmm_multip function to apply Tmat matrix using fmm
% It calculates Tmat*density, where
% zmat = zeros(nnn);
% Tmat = [zmat,        zmat,          D_p2_to_p1,  -S_p2_to_p1  ... D_pn_to_p1,  -S_pn_to_p1;
%         zmat,        zmat,          Dp_p2_to_p1, -Sp_p2_to_p1 ... Dp_pn_to_p1, -Sp_pn_to_p1;
%         D_p1_to_p2,  -S_p1_to_p2,   zmat,        zmat         ... D_pn_to_p2,  -S_pn_to_p2;
%         Dp_p1_to_p2, -Sp_p1_to_p2,  zmat,        zmat         ... Dp_pn_to_p2, -Sp_pn_to_p2
%         ...          ...            ...          ...          ... ...          ...
%         D_p1_to_pn,  -S_p1_to_pn,   D_p2_to_pn,  -S_p2_to_pn  ... zmat,        zmat;
%         Dp_p1_to_pn, -Sp_p1_to_pn,  Dp_p2_to_pn, -Sp_p2_to_pn ... zmat,        zmat];
% 
% where the layer potentials here are the layered medium Green's function.
% For the freespace part this code checkout Tmat_fast_matvec
%
% input : 
%     eps - precision requested
%     zks  - complex number, Helmholtz wave numbers
%     somm_disc - information about sommerfeld discretization
%     npanels - number of panels in plist
%     nnn -  number of points in each pinfo
%     plist - list of panels
%     sqpw   - sqrt(pinfo.wts) % circle points
%     density - vector input to multiply
% output :
%     out - output Tmat*density

dens = density .* repmat(sqpw,2*npanels,1);
pstruct = cat(1, plist{:});

srcuse = [];
srcuse.sources = horzcat(pstruct.r);
dipvec = horzcat(pstruct.n);
srcuse.dipvec = dipvec;
[~, ns] = size(srcuse.sources);
srcuse.charges = zeros(ns, 1);
srcuse.dipstr = zeros(ns, 1);

ishifts = repmat(2*nnn*(0:npanels-1), [nnn,1]);
ishifts = ishifts(:);

indd = repmat(1:nnn, [1,npanels]);
indd = indd(:) + ishifts;

inds = repmat((nnn+1):2*nnn, [1,npanels]);
inds = inds(:) + ishifts;


dipstr = dens(indd);
charges = -dens(inds);
srcuse.dipstr = dipstr.';
srcuse.charges = charges.';
fprintf('starting self correction\n')
corr_vec = complex(2*nnn*npanels,1);
if isa(Tself, 'cell')
    for i=1:npanels
        istart = (i-1)*2*nnn + 1;
        iend = i*2*nnn;
        corr_vec(istart:iend) = Tself{i}*density(istart:iend);
    end
else
    duse = reshape(density, [2*nnn, npanels]);
    corr_vec = Tself*duse;
    corr_vec = reshape(corr_vec, 2*nnn*npanels,1);
end
        
fprintf('starting fmm\n')
fprintf('eps = %d\n', eps);
pg  = 2;
U = hfmm2d(eps, zks(1), srcuse, pg);
sol1 = U.pot(:).*repmat(sqpw, npanels, 1);
sol2 = U.grad(1,:).*srcuse.dipvec(1,:) + U.grad(2,:).*srcuse.dipvec(2,:);
sol2 = sol2(:).*repmat(sqpw,npanels,1);

out = complex(zeros(size(density)));
out(indd) = sol1;
out(inds) = sol2;

% Now compute the correction to the layered medium
ymin = min(srcuse.sources(2,:));
ymax = max(srcuse.sources(2,:));

nleg = 24;
[x, ~, umat, ~] = lege.exps(nleg);
yvals = srcuse.sources(2,:);
yvals = yvals(:);

xvals = srcuse.sources(1,:);
xvals = xvals(:);
yy = (yvals - ymin)/(ymax - ymin)*2 - 1;
pols = lege.pols(yy, nleg-1);
int_mat = pols.'*umat;

yleg  = (x+1)/2*(ymax-ymin) + ymin;

charges_anterp = int_mat.*charges;
dipstr_anterp = int_mat.*dipstr;

pot = complex(zeros(ns,nleg));
gradx = complex(zeros(ns, nleg));
grady = complex(zeros(ns, nleg));

xfac1 = somm_disc.xfac1;
yfac1 = somm_disc.yfac1;
yfac2 = somm_disc.yfac2;
w1 = somm_disc.w1;
fprintf('starting layered medium\n')

rfac = (zks(1)^2 - zks(2)^2)./(1j*yfac1.*(1j*yfac1 + 1j*yfac2).^2).*w1;
for ii = 1:nleg
    ytarg = yleg(ii);
    for jj = 1:nleg
        ysrc = yleg(jj);
        zexp = exp(1j*yfac1*(ysrc + ytarg));
        cuse = charges_anterp(:,jj);
        dxuse = dipstr_anterp(:,jj).*dipvec(1,:).';
        dyuse = dipstr_anterp(:,jj).*dipvec(2,:).';
        f1lm = finufft1d3(xvals, cuse, -1, 1e-13, xfac1);
        f2lm = finufft1d3(xvals, dxuse, -1, 1e-13, xfac1);
        f3lm = finufft1d3(xvals, dyuse, -1, 1e-13, xfac1);
        
        fuse = f1lm -1j*xfac1.*f2lm + 1j*yfac1.*f3lm;
        fuse = fuse.*zexp;
        rf = rfac.*fuse;
        rfx = rfac.*fuse.*(1j*xfac1);
        rfy = rfac.*fuse.*(1j*yfac1);
        pot(:,ii) = pot(:,ii) + finufft1d3(xfac1, rf, 1, 1e-13, xvals);
        gradx(:,ii) = gradx(:,ii) + finufft1d3(xfac1, rfx, 1, 1e-13, xvals);
        grady(:,ii) = grady(:,ii) + finufft1d3(xfac1, rfy, 1, 1e-13, xvals);

    end
end


pot = sum(pot.*int_mat, 2);
gradx = sum(gradx.*int_mat, 2);
grady = sum(grady.*int_mat, 2);
nx = srcuse.dipvec(1,:).';
ny = srcuse.dipvec(2,:).';
gradn = gradx.*nx + grady.*ny;

potc = complex(zeros(size(out)));
potc(indd) = pot.*repmat(sqpw, npanels, 1);
potc(inds) = gradn.*repmat(sqpw, npanels, 1);

out = out(:) + potc(:) - corr_vec(:);

end