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
mean_first = 70001; %50001;
mean_last = 6000000;

%SNAG
MSE_SNAG = zeros(1,runs);
for r=1:runs
    iterates_SNAG = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x + alpha*(x - xm1) - gamma*((X(:,k).'*(x + alpha*(x - xm1)))*X(:,k)-y(k)*X(:,k));
        iterates_SNAG(:,i+1) = x;
        xm1 = iterates_SNAG(:,i);
    end

    iterates_centered_SNAG = iterates_SNAG-xsol;
    iterates_dist_SNAG = vecnorm(iterates_centered_SNAG);
    MSE_SNAG(r) = mean(iterates_dist_SNAG(mean_first:mean_last).^2);
end
mean_MSE_SNAG = mean(MSE_SNAG)
std_MSE_SNAG = std(MSE_SNAG)/sqrt(runs)

%SNAG_approx
MSE_SNAG_approx = zeros(1,runs);
for r=1:runs
    iterates_SNAG_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG_approx(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x + alpha*(x - xm1) - gamma*((X(:,k).'*xsol)*X(:,k)-y(k)*X(:,k) + A*(x + alpha*(x - xm1) - xsol));
        iterates_SNAG_approx(:,i+1) = x;
        xm1 = iterates_SNAG_approx(:,i);
    end

    iterates_centered_SNAG_approx = iterates_SNAG_approx-xsol;
    iterates_dist_SNAG_approx = vecnorm(iterates_centered_SNAG_approx);
    MSE_SNAG_approx(r) = mean(iterates_dist_SNAG_approx(mean_first:mean_last).^2);
end
mean_MSE_SNAG_approx = mean(MSE_SNAG_approx)
std_MSE_SNAG_approx = std(MSE_SNAG_approx)/sqrt(runs)

%SNAG_RR
MSE_SNAG_RR = zeros(1,runs);
for r=1:runs
    iterates_SNAG_RR = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG_RR(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*(x + alpha*(x - xm1)))*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SNAG_RR(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SNAG_RR(:,(i-1)*n+j);
        end
    end
    iterates_centered_SNAG_RR = iterates_SNAG_RR-xsol;
    iterates_dist_SNAG_RR = vecnorm(iterates_centered_SNAG_RR);
    MSE_SNAG_RR(r) = mean(iterates_dist_SNAG_RR(mean_first:mean_last).^2);
end
mean_MSE_SNAG_RR = mean(MSE_SNAG_RR)
std_MSE_SNAG_RR = std(MSE_SNAG_RR)/sqrt(runs)

%SNAG_RR_approx
MSE_SNAG_RR_approx = zeros(1,runs);
for r=1:runs
    iterates_SNAG_RR_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG_RR_approx(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x + alpha*(x - xm1) - xsol));
            iterates_SNAG_RR_approx(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SNAG_RR_approx(:,(i-1)*n+j);
        end
    end
    iterates_centered_SNAG_RR_approx = iterates_SNAG_RR_approx-xsol;
    iterates_dist_SNAG_RR_approx = vecnorm(iterates_centered_SNAG_RR_approx);
    MSE_SNAG_RR_approx(r) = mean(iterates_dist_SNAG_RR_approx(mean_first:mean_last).^2);
end
mean_MSE_SNAG_RR_approx = mean(MSE_SNAG_RR_approx)
std_MSE_SNAG_RR_approx = std(MSE_SNAG_RR_approx)/sqrt(runs)

%SNAG_SO
cycles_per_perm = 150000/n;
MSE_SNAG_SO = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SNAG_SO = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG_SO(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*(x + alpha*(x - xm1)))*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SNAG_SO(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SNAG_SO(:,(i-1)*n+j);
        end
    end
    iterates_centered_SNAG_SO = iterates_SNAG_SO-xsol;
    iterates_dist_SNAG_SO = vecnorm(iterates_centered_SNAG_SO);
    for i=1:(niter/(n*cycles_per_perm))
        MSE_SNAG_SO((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SNAG_SO(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
end
mean_MSE_SNAG_SO = mean(MSE_SNAG_SO)
std_MSE_SNAG_SO = std(MSE_SNAG_SO)/sqrt(runs*niter/(n*cycles_per_perm))

%SNAG_SO_approx
MSE_SNAG_SO_approx = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SNAG_SO_approx = zeros(d,niter+1);
    x = xsol;
    xm1 = x;
    iterates_SNAG_SO_approx(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x + alpha*(x - xm1) - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x + alpha*(x - xm1) - xsol));
            iterates_SNAG_SO_approx(:,(i-1)*n+j+1) = x;
            xm1 = iterates_SNAG_SO_approx(:,(i-1)*n+j);
        end
    end
    iterates_centered_SNAG_SO_approx = iterates_SNAG_SO_approx-xsol;
    iterates_dist_SNAG_SO_approx = vecnorm(iterates_centered_SNAG_SO_approx);
    for i=1:(niter/(n*cycles_per_perm))
        MSE_SNAG_SO_approx((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SNAG_SO_approx(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
end
mean_MSE_SNAG_SO_approx = mean(MSE_SNAG_SO_approx)
std_MSE_SNAG_SO_approx = std(MSE_SNAG_SO_approx)/sqrt(runs*niter/(n*cycles_per_perm))

%variance_A
variance_A = zeros(d);
for i=1:n
    variance_A = variance_A + ((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b))*((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b)).'/n;
end

%variance_SNAG
variance_SNAG_approx = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    variance_SNAG_approx = variance_SNAG_approx + gamma*sigma2/(2*lambda*(1-alpha));
end
disc_points = 5000000;
variance_SNAG = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    xi = (1+alpha)*(1-gamma*lambda);
    chi = alpha*(1-gamma*lambda);
    transfer_sq_SNAG = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_SNAG_i = 0;
    for j=1:disc_points
        variance_SNAG_i = variance_SNAG_i + sigma2*transfer_sq_SNAG(j/disc_points)/disc_points;
    end
    variance_SNAG = variance_SNAG + variance_SNAG_i;
end
variance_SNAG_approx
variance_SNAG
mean_MSE_SNAG_approx

%variance_SNAG_SO
variance_SNAG_SO = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    xi = (1+alpha)*(1-gamma*lambda);
    chi = alpha*(1-gamma*lambda);
    transfer_sq_SNAG = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_SO_i = 0;
    for j=1:(n-1)
        variance_SO_i = variance_SO_i + sigma2*transfer_sq_SNAG(j/n)/(n-1);
    end
    variance_SNAG_SO = variance_SNAG_SO + variance_SO_i;
end
variance_SNAG_SO
mean_MSE_SNAG_SO_approx

%variance_SNAG_RR
variance_SNAG_RR = 0;
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
    xi = (1+alpha)*(1-gamma*lambda);
    chi = alpha*(1-gamma*lambda);
    transfer_sq_SNAG = @(f) (gamma^2/(1+xi^2+chi^2-2*xi*(1+chi)*cos(2*pi*f)+2*chi*cos(4*pi*f)));
    variance_RR_i = 0;
    for j=1:disc_points
        variance_RR_i = variance_RR_i + sigma2*transfer_sq_SNAG(j/disc_points)*spectral_density_RR_2(j/disc_points)/disc_points;
    end
    variance_SNAG_RR = variance_SNAG_RR + variance_RR_i;
end
variance_SNAG_RR
mean_MSE_SNAG_RR_approx

clear iterates_centered_SNAG iterates_centered_SNAG_approx 
clear iterates_centered_SNAG_RR iterates_centered_SNAG_RR_approx
clear iterates_centered_SNAG_SO iterates_centered_SNAG_SO_approx 
clear iterates_SNAG iterates_SNAG_approx 
clear iterates_SNAG_RR iterates_SNAG_RR_approx
clear iterates_SNAG_SO iterates_SNAG_SO_approx