function [u, gradu] = eval_lm_targ_fmm(eps, zks, somm_disc, npanels, nnn, plist, sqpw, targinfo, density)
%
% This subroutine evaluates \sum D_{lm} [\sigma] -S_{lm} [\mu] at a
% collection of targets
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
%     targinfo - target info struct
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
    
pg = 0;
pgt  = 2;
U = hfmm2d(eps, zks(1), srcuse, pg, targinfo.r, pgt);

% Now compute the correction to the layered medium
ymin = min(srcuse.sources(2,:));
ymax = max(srcuse.sources(2,:));

ymintarg = min(targinfo.r(2,:));
ymaxtarg = max(targinfo.r(2,:));

ymin = min(ymin, ymintarg);
ymax = max(ymax, ymaxtarg);

nleg = 24;
[x, ~, umat, ~] = lege.exps(nleg);
yvals = srcuse.sources(2,:);
yvals = yvals(:);

xvals = srcuse.sources(1,:);
xvals = xvals(:);

yvalst = targinfo.r(2,:);
yvalst = yvalst(:);

xvalst = targinfo.r(1,:);
xvalst = xvalst(:);

yy = (yvals - ymin)/(ymax - ymin)*2 - 1;
yyt = (yvalst - ymin)/(ymax - ymin)*2 - 1;

pols = lege.pols(yy, nleg-1);
polst = lege.pols(yyt, nleg-1);

int_mat = pols.'*umat;
int_matt = polst.'*umat;

yleg  = (x+1)/2*(ymax-ymin) + ymin;

charges_anterp = int_mat.*charges;
dipstr_anterp = int_mat.*dipstr;

[~, nt] = size(targinfo.r);

pot = complex(zeros(nt, nleg));
gradx = complex(zeros(nt, nleg));
grady = complex(zeros(nt, nleg));

xfac1 = somm_disc.xfac1;
yfac1 = somm_disc.yfac1;
yfac2 = somm_disc.yfac2;
w1 = somm_disc.w1;


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
        pot(:,ii) = pot(:,ii) + finufft1d3(xfac1, rf, 1, 1e-13, xvalst);
        gradx(:,ii) = gradx(:,ii) + finufft1d3(xfac1, rfx, 1, 1e-13, xvalst);
        grady(:,ii) = grady(:,ii) + finufft1d3(xfac1, rfy, 1, 1e-13, xvalst);

    end
end


pot = sum(pot.*int_matt, 2);
gradx = sum(gradx.*int_matt, 2);
grady = sum(grady.*int_matt, 2);

pot1 = U.pottarg(:);
gradx1 = U.gradtarg(1,:); gradx1 = gradx1(:);
grady1 = U.gradtarg(2,:); grady1 = grady1(:);

u = pot + pot1;
gradu = complex(zeros(2,nt));
gradu(1,:) = gradx + gradx1;
gradu(2,:) = grady + grady1;

end