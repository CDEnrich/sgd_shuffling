figure('position',[1,1,1200,350],'DefaultAxesFontSize',14);
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(iterates_dist_SGD(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGD_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGD with replacement','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SGD_RR(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGD_RR_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGD-RR','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SGD_SO(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGD_SO_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGD-SO','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)
ylim([0 2.5e-08])

figure('position',[1,1,1200,350],'DefaultAxesFontSize',14);
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(iterates_dist_SGDM(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGDM_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGDM with replacement','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SGDM_RR(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGDM_RR_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGDM-RR','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SGDM_SO(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SGDM_SO_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SGDM-SO','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)
ylim([0 1e-07])

figure('position',[1,1,1200,350],'DefaultAxesFontSize',14);
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(iterates_dist_SNAG(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SNAG_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SNAG with replacement','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SNAG_RR(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SNAG_RR_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SNAG-RR','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)

nexttile
plot(iterates_dist_SNAG_SO(mean_first:mean_last).^2)
hold on
plot(iterates_dist_SNAG_SO_approx(mean_first:mean_last).^2)
legend('Original','0-th order model','FontSize',17)
title('SNAG-SO','FontSize',17)
xlabel('Iteration','FontSize',17)
ylabel('Squared distance to optimum','FontSize',17)
ylim([0 1e-07])