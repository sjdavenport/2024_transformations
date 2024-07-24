FWHM_vec = 0:5;
fpr_orig = zeros(1, length(FWHM_vec));
fpr_apower = zeros(1, length(FWHM_vec));
fpr_asinh_clt = zeros(1, length(FWHM_vec));
fpr_asinh_exact = zeros(1, length(FWHM_vec));
saveloc = '/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_transformations/';
% fpr_asinh_exact = 0;

nsim = 5000;

nsubj = 40;
for J = 1:length(FWHM_vec)
    for I = 1:nsim
        I
        MNImask = imgload('MNImask') > 0;
        slice = 50;
        MNImask2D = MNImask(:,:,slice);
        noise = trnd(3, [size(MNImask2D), nsubj]);
        smooth_noise = fast_conv(noise, FWHM, 2);
        data = vec_data( smooth_noise, MNImask2D );

        data_transformed_apower = apower(data);
        data_transformed_asinh = asinh(data./std(data,0,2));

        tstat_orig = mvtstat(data);
        tstat_apower = mvtstat(data_transformed_apower);
        tstat_asinh = mvtstat(data_transformed_asinh);

        orig_thresh = fastperm( data, nsubj, alpha, nboot, show_loader);
        apower_thresh = fastperm( data_transformed_apower, nsubj, alpha, nboot, show_loader);
        asinh_thresh_clt = fastperm( data_transformed_asinh, nsubj, alpha, nboot, show_loader);
        % asinh_thresh_exact = perm_asinh( data, [sum(MNImask2D(:)),1] > 0, alpha, nboot, show_loader);

        fpr_orig(J) = fpr_orig(J) + (max(tstat_orig) > orig_thresh);
        fpr_apower(J) = fpr_apower(J) + (max(tstat_apower) > apower_thresh);
        fpr_asinh_clt(J) = fpr_asinh_clt(J) + (max(tstat_asinh) > asinh_thresh_clt);

        asinh_thresh_exact = perm_asinh( data, alpha, nboot, show_loader);
        fpr_asinh_exact = fpr_asinh_exact + (max(tstat_asinh) > asinh_thresh_exact);
        % fpr_asinh_exact = fpr_asinh_exact + (max(tstat_asinh) > asinh_thresh_exact);
    end
    save([saveloc, 'nsubj_', num2str(nsubj)], 'fpr_orig', 'fpr_apower', 'fpr_asinh_clt', 'fpr_asinh_exact')
end
fpr_orig = fpr_orig/nsim;
fpr_apower = fpr_apower/nsim;
fpr_asinh_clt = fpr_asinh_clt/nsim;
fpr_asinh_exact = fpr_asinh_exact/nsim;
save([saveloc, 'nsubj_', num2str(nsubj)], 'fpr_orig', 'fpr_apower', 'fpr_asinh_clt', 'fpr_asinh_exact')

%%
fpr_orig/nsim
fpr_apower/nsim
fpr_asinh_clt/nsim
fpr_asinh_exact/nsim

%%
fpr_asinh_exact = 0;

nsim = 1000;

alpha = 0.05;

FWHM = 4;
nsubj = 30;
for I = 1:nsim
    I
    MNImask = imgload('MNImask') > 0;
    slice = 50;
    MNImask2D = MNImask(:,:,slice);
    noise = trnd(3, [size(MNImask2D), nsubj]);
    smooth_noise = fast_conv(noise, FWHM, 2);
    data = vec_data( smooth_noise, MNImask2D );

    data_transformed_asinh = asinh(data./std(data,0,2));

    tstat_asinh = mvtstat(data_transformed_asinh);

    asinh_thresh_exact = perm_asinh( data, alpha, nboot, show_loader);

    fpr_asinh_exact = fpr_asinh_exact + (max(tstat_asinh) > asinh_thresh_exact);
end

fpr_asinh_exact/nsim