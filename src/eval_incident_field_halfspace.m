function [u, gradu, splits] = eval_incident_field_halfspace(zks, somm_disc, src_info, targ, strengths, type)

if nargin < 6
    type = 's';
end

if nargin < 5
    [~, ns] = size(src_info.r);
    strengths = ones(ns,1);
end

s_info = [];
s_info.r = src_info.r;
sk1 = kernel('helm', 's', zks(1));
[~, nt] = size(targ.r);

if strcmpi(type, 's')

    src = src_info.r;
    trg = targ.r;
    [u1, grad] = chnk.helm2d.green(zks(1), src, trg);
    u1 = u1*strengths(:);
    gx1 = grad(:,:,1)*strengths(:);
    gy1 = grad(:,:,2)*strengths(:);
    
    s_info.charges = 1;
    [ucorr, graducorr] = eval_sommerfeld_correction(zks, somm_disc, s_info, targ);
    u = u1 + ucorr*strengths(:);
    if nargout > 1
        gxcorr = graducorr(:,:,1)*strengths(:);
        gycorr = graducorr(:,:,2)*strengths(:);
        gradu = complex(zeros(2,nt));
        gradu(1,:) = gx1 + gxcorr;
        gradu(2,:) = gy1 + gycorr;        
    end
end

if nargout > 2
    splits = [];
    splits.u = u1;
    splits.gradu = complex(zeros(2,nt));
    splits.gradu(1,:) = gx1;
    splits.gradu(2,:) = gy1;
    splits.ucorr = ucorr*strengths(:);
    
    splits.graducorr = complex(zeros(2,nt));
    splits.graducorr(1,:) = gxcorr;
    splits.graducorr(2,:) = gycorr;
end


end