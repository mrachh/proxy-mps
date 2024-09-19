function [u, gradu, mats] = eval_sommerfeld_correction(zks, somm_disc, src_info, targ_info)
%
%  This subroutine evaluates the correction to free space Green's function
%  in order to compute the layered medium Green's function
%  
%  Currently this subroutine only supports sources and targets in 
%  the upper half plane
%
xfac1 = somm_disc.xfac1;
yfac1 = somm_disc.yfac1;
yfac2 = somm_disc.yfac2;
w1 = somm_disc.w1;

rfac = (zks(1)^2 - zks(2)^2)./(1j*yfac1.*(1j*yfac1 + 1j*yfac2).^2);

sr = src_info.r;
[~, ns] = size(sr);
charges = zeros(1,ns);
dipstr = zeros(1,ns);
dipvec = zeros(2,ns);

if isfield(src_info, 'charges')
    charges = src_info.charges(:).';
end

if isfield(src_info, 'dipstr')
    dipstr = src_info.dipstr(:).';
end

if isfield(src_info, 'dipvec')
    dipvec = src_info.dipvec;
end

dipvec = dipstr.*dipvec;

tr = targ_info.r;
sexp0 = rfac.*exp(1j*(-xfac1*sr(1,:) + yfac1*sr(2,:)));
sexp = sexp0.*charges;
sexp = sexp + sexp0.*(-1j*xfac1*dipvec(1,:) + 1j*yfac1*dipvec(2,:));
texp = (w1.').*exp(1j*(xfac1*tr(1,:) + yfac1*tr(2,:))).';

u = texp*sexp;



if nargout >1
    [nt, ns] = size(u);
    gradu = complex(zeros(nt,ns,2));
    texpx = texp.*1j.*xfac1.';
    texpy = texp.*1j.*yfac1.';
    gradu(:,:,1) = texpx*sexp;
    gradu(:,:,2) = texpy*sexp;  
end

if nargout > 2
    mats = [];
    mats.sexp = sexp;
    mats.texp = texp;
    mats.texpx = texpx;
    mats.texpy = texpy;
end


end