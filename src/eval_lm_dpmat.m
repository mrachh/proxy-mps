function dpmat = eval_lm_dpmat(zks, somm_disc, s, t)
    dkp1 = kernel('helm', 'dp', zks(1));
    dkpmat = dkp1.eval(s, t);
    dkpmat(isnan(dkpmat)) = 0;

    src_info = [];
    src_info.r = s.r(:,:);
    src_info.dipvec = s.n(:,:);
    src_info.dipstr = 1;

    [~, gradscorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
    dcorr = gradscorr(:,:,1).*t.n(1,:).' + gradscorr(:,:,2).*t.n(2,:).';
    dpmat = dkpmat + dcorr;

end