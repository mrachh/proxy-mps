function [amat] = eval_lm_mat(zks, somm_disc, s, t, type, iquad)
    if nargin < 6
        iquad = 0;
    end
    w = s.wts(:);
    src_info = [];
    src_info.r = s.r(:,:);
    src_info.dipvec = s.n(:,:);


    if strcmpi(type, 's')
        K = kernel('helm', 's', zks(1));
        if iquad == 1
            amat0 = chunkermat(s, K);
        else
            amat0 = K.eval(s,t).*(w.');
        end
        src_info.charges = 1;

        acorr = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
        amat = amat0 + acorr.*(w.');
        

    elseif strcmpi(type, 'd')
        K = kernel('helm', 'd', zks(1));
        if iquad == 1
            amat0 = chunkermat(s, K);
        else
            amat0 = K.eval(s,t).*(w.');
        end
        src_info.dipstr = 1;

        [acorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
        amat = amat0 + acorr.*(w.');

    elseif strcmpi(type, 'sp')
        K = kernel('helm', 'sp', zks(1));
        if iquad == 1
            amat0 = chunkermat(s, K);
        else
            amat0 = K.eval(s,t).*(w.');
        end
        src_info.charges = 1;

        [~, gradscorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
        acorr = gradscorr(:,:,1).*(t.n(1,:).') + gradscorr(:,:,2).*(t.n(2,:).');
        amat = amat0 + acorr.*(w.');


    elseif strcmpi(type, 'dp')
        K = kernel('helm', 'dp', zks(1));
        if iquad == 1
            amat0 = chunkermat(s, K);
        else
            amat0 = K.eval(s,t).*(w.');
        end
        src_info.dipstr = 1;

        [~, gradscorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
        acorr = gradscorr(:,:,1).*(t.n(1,:).') + gradscorr(:,:,2).*(t.n(2,:).');
        amat = amat0 + acorr.*(w.');
    end

end