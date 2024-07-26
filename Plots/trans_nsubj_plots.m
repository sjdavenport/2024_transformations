a = load('/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_transformations/nsubj_results/nvox_1000.mat');

nsubj_vec = 10:10:100;

alpha = 0.05;
nsim = 5000;



plot(nsubj_vec, a.fpr_apower(:,1)/5000, 'LineWidth', 6 );
hold on
plot(nsubj_vec, a.fpr_asinh_clt(:,1)/5000, 'LineWidth', 6 );
plot(nsubj_vec, a.fpr_asinh_exact(:,1)/5000, 'LineWidth', 6 );
% plot(nsubj_vec, a.fpr_orig(:,1)/5000, 'LineWidth', 6 );

yline(alpha, '-', 'LineWidth', 3 );
interval = bernstd( alpha, nsim );
yline(interval(1), '--', 'LineWidth', 3 );
yline(interval(2), '--', 'LineWidth', 3 );

ylim([0,0.12])
xlim([10,100])

legend('Apower', 'Asinh clt', 'Asinh exact')

xlabel('Number of subjects')
ylabel('Coverage Rate')

BigFont(30)

matniceplot
fullscreen
saveim('fpr')