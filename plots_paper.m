plot_iterations = 599999;

f1 = figure('position',[1,1,1200,350])
subplot(1,3,1)
plot(iterates_dist_SGD(mean_last-plot_iterations:mean_last).^2) 
hold on
plot(iterates_dist_SGD_approx(mean_last-plot_iterations:mean_last).^2)
legend('Original', '0-th order model')
title('SGD with replacement')
xlabel('Iteration')
ylabel('Squared distance to optimum')

subplot(1,3,2)
plot(iterates_dist_SGD_RR(mean_last-plot_iterations:mean_last).^2) 
hold on
plot(iterates_dist_SGD_RR_approx(mean_last-plot_iterations:mean_last).^2)
legend('Original', '0-th order model')
title('SGD with replacement')
xlabel('Iteration')
ylabel('Squared distance to optimum')

subplot(1,3,3)
sqiterates_dist_SGD_SO = downsample(iterates_dist_SGD_SO(mean_first:mean_last).^2,10);
sqiterates_dist_SGD_SO_approx = downsample(iterates_dist_SGD_SO_approx(mean_first:mean_last).^2,10);
plot(sqiterates_dist_SGD_SO) 
hold on
plot(sqiterates_dist_SGD_SO_approx)
legend('Original', '0-th order model')
title('SGD with replacement')
xlabel('Downsampled iteration')
ylabel('Squared distance to optimum')