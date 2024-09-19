function out = Tmat_fast_matvec(eps, zk, npanels, nnn, plist, sqpw, Tself, density)
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
srcuse.dipvec = horzcat(pstruct.n);
[~, ns] = size(srcuse.sources);
srcuse.charges = zeros(ns, 1);
srcuse.dipstr = zeros(ns, 1);

ishifts = repmat(2*nnn*(0:npanels-1), [nnn,1]);
ishifts = ishifts(:);

indd = repmat(1:nnn, [1,npanels]);
indd = indd(:) + ishifts;

inds = repmat((nnn+1):2*nnn, [1,npanels]);
inds = inds(:) + ishifts;

srcuse.dipstr = dens(indd).';
srcuse.charges = -dens(inds).';


duse = reshape(density, [2*nnn, npanels]);
corr_vec = Tself*duse;
corr_vec = reshape(corr_vec, 2*nnn*npanels,1);
        

pg  = 2;
U = hfmm2d(eps, zk, srcuse, pg);
sol1 = U.pot(:).*repmat(sqpw, npanels, 1);
sol2 = U.grad(1,:).*srcuse.dipvec(1,:) + U.grad(2,:).*srcuse.dipvec(2,:);
sol2 = sol2(:).*repmat(sqpw,npanels,1);

out = complex(zeros(size(density)));
out(indd) = sol1;
out(inds) = sol2;

out = out(:) - corr_vec(:);

return