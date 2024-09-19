function dmat = eval_lm_dmat(zks, somm_disc, s, t)
    dk1 = kernel('helm', 'd', zks(1));
    dkmat = dk1.eval(s, t);
    dkmat(isnan(dkmat)) = 0;

    src_info = [];
    src_info.r = s.r(:,:);
    src_info.dipvec = s.n(:,:);
    src_info.dipstr = 1;

    [dcorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
    dmat = dkmat + dcorr;

end