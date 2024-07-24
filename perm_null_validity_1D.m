nsubj_vec = 10:10:100;

fpr_orig = zeros(length(nsubj_vec), 3);
fpr_apower = zeros(length(nsubj_vec), 3);
fpr_asinh_clt = zeros(length(nsubj_vec), 3);
fpr_asinh_exact = zeros(length(nsubj_vec), 3);
saveloc = '/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_transformations/nsubj_results/';
% fpr_asinh_exact = 0;

nsim = 5000;
alpha = 0.05;
nboot = 1000;
show_loader = 0;
nvox = 1000;
for K = 1:3
    for J = 1:length(nsubj_vec)
        nsubj = nsubj_vec(J);
        for I = 1:nsim
            I
            if K == 1
                data = trnd(3, [nvox, nsubj]);
            elseif K == 2
                data = rlap( 3, [nvox, nsubj] );
            elseif K == 3
                data = normrnd(0,1,[nvox, nsubj]);
            end

            data_transformed_apower = apower(data);
            data_transformed_asinh = asinh(data./std(data,0,2));

            tstat_orig = mvtstat(data);
            tstat_apower = mvtstat(data_transformed_apower);
            tstat_asinh = mvtstat(data_transformed_asinh);

            orig_thresh = fastperm( data, nsubj, alpha, nboot, show_loader);
            apower_thresh = fastperm( data_transformed_apower, nsubj, alpha, nboot, show_loader);
            asinh_thresh_clt = fastperm( data_transformed_asinh, nsubj, alpha, nboot, show_loader);
            % asinh_thresh_exact = perm_asinh( data, [sum(MNImask2D(:)),1] > 0, alpha, nboot, show_loader);

            fpr_orig(J, K) = fpr_orig(J, K) + (max(tstat_orig) > orig_thresh);
            fpr_apower(J, K) = fpr_apower(J,K) + (max(tstat_apower) > apower_thresh);
            fpr_asinh_clt(J,K) = fpr_asinh_clt(J,K) + (max(tstat_asinh) > asinh_thresh_clt);

            asinh_thresh_exact = perm_asinh( data, alpha, nboot, show_loader);
            fpr_asinh_exact(J,K) = fpr_asinh_exact(J,K) + (max(tstat_asinh) > asinh_thresh_exact);

            save([saveloc, 'nvox_', num2str(nvox)], 'fpr_orig', 'fpr_apower', 'fpr_asinh_clt', 'fpr_asinh_exact')
        end
    end
end
fpr_orig = fpr_orig/nsim
fpr_apower = fpr_apower/nsim
fpr_asinh_clt = fpr_asinh_clt/nsim
fpr_asinh_exact = fpr_asinh_exact/nsim
save([saveloc, 'nvox_', num2str(nvox)], 'fpr_orig', 'fpr_apower', 'fpr_asinh_clt', 'fpr_asinh_exact')

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