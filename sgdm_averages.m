n = 10; %1000;
d = 5;
gamma = 0.002; %0.0005; %0.01;
alpha = 0.8; %0.9;
runs = 10;

rng(38)

%matrix X and A
X = normrnd(0,1/sqrt(d),[d,n]);
A = X*X.'/n;
eigenvalues = eig(A);
[eigenvectors,D] = eig(A);

%target
e1 = zeros(1,d);
e1(1) = 1;
y = e1*X + normrnd(0,0.1,[1,n]);
b = -sum(y.*X,2)/n;
xsol = -A\b;

niter = 6000000;
mean_first = 120001; %50001;
mean_last = 6000000;

%SGDM
MSE_SGDM = zeros(1,runs);
for r=1:runs
    iterates_SGDM = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x + alpha*(x - xm1) - gamma*((X(:,k).'*x)*X(:,k)-y(k)*X(:,k));
        iterates_SGDM(:,i+1) = x;
        xm1 = iterates_SGDM(:,i);
    end

    iterates_centered_SGDM = iterates_SGDM-xsol;
    iterates_dist_SGDM = vecnorm(iterates_centered_SGDM);
    MSE_SGDM(r) = mean(iterates_dist_SGDM(mean_first:mean_last).^2);
end
mean_MSE_SGDM = mean(MSE_SGDM)
std_MSE_SGDM = std(MSE_SGDM)/sqrt(runs)

%SGDM_approx
MSE_SGDM_approx = zeros(1,runs);
for r=1:runs
    iterates_SGDM_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM_approx(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x + alpha*(x - xm1) - gamma*((X(:,k).'*xsol)*X(:,k)-y(k)*X(:,k) + A*(x - xsol));
        iterates_SGDM_approx(:,i+1) = x;
        xm1 = iterates_SGDM_approx(:,i);
    end

    iterates_centered_SGDM_approx = iterates_SGDM_approx-xsol;
    iterates_dist_SGDM_approx = vecnorm(iterates_centered_SGDM_approx);
    MSE_SGDM_approx(r) = mean(iterates_dist_SGDM_approx(mean_first:mean_last).^2);
end
mean_MSE_SGDM_approx = mean(MSE_SGDM_approx)
std_MSE_SGDM_approx = std(MSE_SGDM_approx)/sqrt(runs)

%SGDM_RR
MSE_SGDM_RR = zeros(1,runs);
for r=1:runs
    iterates_SGDM_RR = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM_RR(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*x)*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SGDM_RR(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SGDM_RR(:,(i-1)*n+j);
        end
    end
    iterates_centered_SGDM_RR = iterates_SGDM_RR-xsol;
    iterates_dist_SGDM_RR = vecnorm(iterates_centered_SGDM_RR);
    MSE_SGDM_RR(r) = mean(iterates_dist_SGDM_RR(mean_first:mean_last).^2);
end
mean_MSE_SGDM_RR = mean(MSE_SGDM_RR)
std_MSE_SGDM_RR = std(MSE_SGDM_RR)/sqrt(runs)

%SGDM_RR_approx
MSE_SGDM_RR_approx = zeros(1,runs);
for r=1:runs
    iterates_SGDM_RR_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM_RR_approx(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x - xsol));
            iterates_SGDM_RR_approx(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SGDM_RR_approx(:,(i-1)*n+j);
        end
    end
    iterates_centered_SGDM_RR_approx = iterates_SGDM_RR_approx-xsol;
    iterates_dist_SGDM_RR_approx = vecnorm(iterates_centered_SGDM_RR_approx);
    MSE_SGDM_RR_approx(r) = mean(iterates_dist_SGDM_RR_approx(mean_first:mean_last).^2);
end
mean_MSE_SGDM_RR_approx = mean(MSE_SGDM_RR_approx)
std_MSE_SGDM_RR_approx = std(MSE_SGDM_RR_approx)/sqrt(runs)

%SGDM_SO
cycles_per_perm = 150000/n;
MSE_SGDM_SO = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SGDM_SO = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM_SO(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*x)*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SGDM_SO(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SGDM_SO(:,(i-1)*n+j);
        end
    end
    iterates_centered_SGDM_SO = iterates_SGDM_SO-xsol;
    iterates_dist_SGDM_SO = vecnorm(iterates_centered_SGDM_SO);
    for i=1:(niter/(n*cycles_per_perm))
        MSE_SGDM_SO((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SGDM_SO(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
end
mean_MSE_SGDM_SO = mean(MSE_SGDM_SO)
std_MSE_SGDM_SO = std(MSE_SGDM_SO)/sqrt(runs*niter/(n*cycles_per_perm))

%SGDM_SO_approx
MSE_SGDM_SO_approx = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SGDM_SO_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SGDM_SO_approx(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x - xsol));
            iterates_SGDM_SO_approx(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SGDM_SO_approx(:,(i-1)*n+j);
        end
    end
    iterates_centered_SGDM_SO_approx = iterates_SGDM_SO_approx-xsol;
    iterates_dist_SGDM_SO_approx = vecnorm(iterates_centered_SGDM_SO_approx);
    for i=1:(niter/(n*cycles_per_perm))
        MSE_SGDM_SO_approx((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SGDM_SO_approx(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
end
mean_MSE_SGDM_SO_approx = mean(MSE_SGDM_SO_approx)
std_MSE_SGDM_SO_approx = std(MSE_SGDM_SO_approx)/sqrt(runs*niter/(n*cycles_per_perm))

%variance_A
variance_A = zeros(d);
for i=1:n
    variance_A = variance_A + ((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b))*((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b)).'/n;
end

%variance_SGDM
variance_SGDM_approx = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    variance_SGDM_approx = variance_SGDM_approx + gamma*sigma2/(2*lambda*(1-alpha));
end
disc_points = 5000000;
variance_SGDM = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    xi = 1+alpha-gamma*lambda;
    chi = alpha;
    transfer_sq_SGDM = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_SGDM_i = 0;
    for j=1:disc_points
        variance_SGDM_i = variance_SGDM_i + sigma2*transfer_sq_SGDM(j/disc_points)/disc_points;
    end
    variance_SGDM = variance_SGDM + variance_SGDM_i;
end
variance_SGDM_approx
variance_SGDM
mean_MSE_SGDM_approx

%variance_SGDM_SO
variance_SGDM_SO = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    xi = 1+alpha-gamma*lambda;
    chi = alpha;
    transfer_sq_SGDM = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_SO_i = 0;
    for j=1:(n-1)
        variance_SO_i = variance_SO_i + sigma2*transfer_sq_SGDM(j/n)/(n-1);
    end
    variance_SGDM_SO = variance_SGDM_SO + variance_SO_i;
end
variance_SGDM_SO
mean_MSE_SGDM_SO_approx

%variance_SGDM_RR
variance_SGDM_RR = 0;
disc_points = 5000000;
cos_sin_term_1 = @(f) ((sin(2*pi*(n-1/2)*f)-sin(pi*(n-1)*f)*cos(pi*n*f))/sin(pi*f));
cos_sin_term_2 = @(f) (sin(pi*(n-1)*f)*sin(pi*n*f)*cos(pi*f)/sin(pi*f)^2);
spectral_density_RR = @(f) (n/(n-1)-(cos_sin_term_1(f)+cos_sin_term_2(f))/(n*(n-1)));
cos_vec = @(f) (cos(2*pi*f*(1:(n-1))));
sum_cos_vec = @(f) (n-(1:(n-1)))*cos_vec(f).';
spectral_density_RR_2 = @(f) (1-2*sum_cos_vec(f)/(n*(n-1)));
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    xi = 1+alpha-gamma*lambda;
    chi = alpha;
    transfer_sq_SGDM = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_RR_i = 0;
    for j=1:disc_points
        variance_RR_i = variance_RR_i + sigma2*transfer_sq_SGDM(j/disc_points)*spectral_density_RR_2(j/disc_points)/disc_points;
    end
    variance_SGDM_RR = variance_SGDM_RR + variance_RR_i;
end
variance_SGDM_RR
mean_MSE_SGDM_RR_approx

clear iterates_centered_SGDM iterates_centered_SGDM_approx 
clear iterates_centered_SGDM_RR iterates_centered_SGDM_RR_approx
clear iterates_centered_SGDM_SO iterates_centered_SGDM_SO_approx 
clear iterates_SGDM iterates_SGDM_approx 
clear iterates_SGDM_RR iterates_SGDM_RR_approx
clear iterates_SGDM_SO iterates_SGDM_SO_approx