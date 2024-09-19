function smat = eval_lm_smat(zks, somm_disc, s, t)
    sk1 = kernel('helm', 's', zks(1));
    skmat = sk1.eval(s, t);
    skmat(isnan(skmat)) = 0;

    src_info = [];
    src_info.r = s.r(:,:);
    src_info.charges = 1;

    [scorr] = eval_sommerfeld_correction(zks, somm_disc, src_info, t);
    smat = skmat + scorr;

end