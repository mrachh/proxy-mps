function spmat = eval_lm_spmat(zks, somm_disc, s, t)
    skp1 = kernel('helm', 'sp', zks(1));
    skpmat = skp1.eval(s, t);
    skpmat(isnan(skpmat)) = 0;

    src_info = [];
    src_info.r = s.r(:,:);
    src_info.charges = 1;

    [~, gradscorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
    scorr = gradscorr(:,:,1).*(t.n(1,:).') + gradscorr(:,:,2).*(t.n(2,:).');
    spmat = skpmat + scorr;

end