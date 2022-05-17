n = 10; %1000;
d = 5;
gamma = 0.002; %0.0005;
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

niter = 6000000; %12000000;
mean_first = 120001; %50001; 
mean_last = 6000000; %12000000; 

%SGD
MSE_SGD = zeros(1,runs);
for r=1:runs
    iterates_SGD = zeros(d,niter+1);
    x = xsol;
    iterates_SGD(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x - gamma*((X(:,k).'*x)*X(:,k)-y(k)*X(:,k));
        iterates_SGD(:,i+1) = x;
    end

    iterates_centered_SGD = iterates_SGD-xsol;
    iterates_dist_SGD = vecnorm(iterates_centered_SGD);
    MSE_SGD(r) = mean(iterates_dist_SGD(mean_first:mean_last).^2);
end
mean_MSE_SGD = mean(MSE_SGD)
std_MSE_SGD = std(MSE_SGD)/sqrt(runs)

%SGD_approx
MSE_SGD_approx = zeros(1,runs);
for r=1:runs
    iterates_SGD_approx = zeros(d,niter+1);
    x = xsol;
    iterates_SGD_approx(:,1) = x;
    for i=1:niter
        k = randi([1,n],1);
        x = x - gamma*((X(:,k).'*xsol)*X(:,k)-y(k)*X(:,k) + A*(x - xsol));
        iterates_SGD_approx(:,i+1) = x;
    end

    iterates_centered_SGD_approx = iterates_SGD_approx-xsol;
    iterates_dist_SGD_approx = vecnorm(iterates_centered_SGD_approx);
    MSE_SGD_approx(r) = mean(iterates_dist_SGD_approx(mean_first:mean_last).^2);
end
mean_MSE_SGD_approx = mean(MSE_SGD_approx)
std_MSE_SGD_approx = std(MSE_SGD_approx)/sqrt(runs)

%SGD_RR
MSE_SGD_RR = zeros(1,runs);
for r=1:runs
    iterates_SGD_RR = zeros(d,niter+1);
    x = xsol;
    iterates_SGD_RR(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x - gamma*((X(:,p(j)).'*x)*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SGD_RR(:,(i-1)*n+j+1) = x;
        end
    end
    
    iterates_centered_SGD_RR = iterates_SGD_RR-xsol;
    iterates_dist_SGD_RR = vecnorm(iterates_centered_SGD_RR);
    MSE_SGD_RR(r) = mean(iterates_dist_SGD_RR(mean_first:mean_last).^2);
end
mean_MSE_SGD_RR = mean(MSE_SGD_RR)
std_MSE_SGD_RR = std(MSE_SGD_RR)/sqrt(runs)

%SGD_RR_approx
MSE_SGD_RR_approx = zeros(1,runs);
for r=1:runs
    iterates_SGD_RR_approx = zeros(d,niter+1);
    x = xsol;
    iterates_SGD_RR_approx(:,1) = x;
    for i=1:(niter/n)
        p = randperm(n);
        for j=1:n
            x = x - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x - xsol));
            iterates_SGD_RR_approx(:,(i-1)*n+j+1) = x;
        end
    end
    iterates_centered_SGD_RR_approx = iterates_SGD_RR_approx-xsol;
    iterates_dist_SGD_RR_approx = vecnorm(iterates_centered_SGD_RR_approx);
    MSE_SGD_RR_approx(r) = mean(iterates_dist_SGD_RR_approx(mean_first:mean_last).^2);
end
mean_MSE_SGD_RR_approx = mean(MSE_SGD_RR_approx)
std_MSE_SGD_RR_approx = std(MSE_SGD_RR_approx)/sqrt(runs)

%SGD_SO
%MSE_SGD_SO = zeros(1,runs);
cycles_per_perm = 150000/n;
MSE_SGD_SO = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SGD_SO = zeros(d,niter+1); 
    x = xsol;
    iterates_SGD_SO(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x - gamma*((X(:,p(j)).'*x)*X(:,p(j))-y(p(j))*X(:,p(j)));
            iterates_SGD_SO(:,(i-1)*n+j+1) = x;
        end
    end
    iterates_centered_SGD_SO = iterates_SGD_SO-xsol;
    iterates_dist_SGD_SO = vecnorm(iterates_centered_SGD_SO);
    %MSE_SGD_SO_r = 0;
    for i=1:(niter/(n*cycles_per_perm))
        %MSE_SGD_SO_r = MSE_SGD_SO_r + mean(iterates_dist_SGD_SO(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
        MSE_SGD_SO((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SGD_SO(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
    %MSE_SGD_SO(r) = MSE_SGD_SO_r/(niter/(n*cycles_per_perm));
end
mean_MSE_SGD_SO = mean(MSE_SGD_SO)
%std_MSE_SGD_SO = std(MSE_SGD_SO)/sqrt(runs)
std_MSE_SGD_SO = std(MSE_SGD_SO)/sqrt(runs*niter/(n*cycles_per_perm))

%SGD_SO_approx
%MSE_SGD_SO_approx = zeros(1,runs);
MSE_SGD_SO_approx = zeros(1,runs*niter/(n*cycles_per_perm));
for r=1:runs
    iterates_SGD_SO_approx = zeros(d,niter+1);
    x = xsol;
    iterates_SGD_SO_approx(:,1) = x;
    p = randperm(n);
    for i=1:(niter/n)
        if mod(i,cycles_per_perm) == 1
            p = randperm(n);
        end
        for j=1:n
            x = x - gamma*((X(:,p(j)).'*xsol)*X(:,p(j))-y(p(j))*X(:,p(j)) + A*(x - xsol));
            iterates_SGD_SO_approx(:,(i-1)*n+j+1) = x;
        end
    end
    iterates_centered_SGD_SO_approx = iterates_SGD_SO_approx-xsol;
    iterates_dist_SGD_SO_approx = vecnorm(iterates_centered_SGD_SO_approx);
    %MSE_SGD_SO_approx_r = 0;
    for i=1:(niter/(n*cycles_per_perm))
        %MSE_SGD_SO_approx_r = MSE_SGD_SO_approx_r + mean(iterates_dist_SGD_SO_approx(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
        MSE_SGD_SO_approx((r-1)*niter/(n*cycles_per_perm)+i) = mean(iterates_dist_SGD_SO_approx(((i-1)*n*cycles_per_perm + mean_first):(i*n*cycles_per_perm)).^2);
    end
    %MSE_SGD_SO_approx(r) = MSE_SGD_SO_approx_r/(niter/(n*cycles_per_perm));
end
mean_MSE_SGD_SO_approx = mean(MSE_SGD_SO_approx)
%std_MSE_SGD_SO_approx = std(MSE_SGD_SO_approx)/sqrt(runs)
std_MSE_SGD_SO_approx = std(MSE_SGD_SO_approx)/sqrt(runs*niter/(n*cycles_per_perm))

%variance_A
variance_A = zeros(d);
for i=1:n
    variance_A = variance_A + ((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b))*((X(:,i).'*xsol)*X(:,i)-y(i)*X(:,i)-(A*xsol+b)).'/n;
end

%variance_SGD
variance_SGD = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    variance_SGD = variance_SGD + gamma*sigma2/(2*lambda - gamma*lambda.^2);
end    
variance_SGD
mean_MSE_SGD_approx

%variance_SGD_SO
variance_SGD_SO = 0;
for i=1:d
    lambda = D(i,i);
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec;
    transfer_sq_SGD = @(f) (gamma^2/(1+(1-gamma*lambda)^2-2*(1-gamma*lambda)*cos(2*pi*f)));
    variance_SO_i = 0;
    for j=1:(n-1)
        variance_SO_i = variance_SO_i + sigma2*transfer_sq_SGD(j/n)/(n-1);
    end
    variance_SGD_SO = variance_SGD_SO + variance_SO_i;
end
variance_SGD_SO
mean_MSE_SGD_SO_approx

%variance_SGD_RR
variance_SGD_RR = 0;
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
    transfer_sq_SGD = @(f) (gamma^2/(1+(1-gamma*lambda)^2-2*(1-gamma*lambda)*cos(2*pi*f)));
    variance_RR_i = 0;
    for j=1:disc_points
        variance_RR_i = variance_RR_i + sigma2*transfer_sq_SGD(j/disc_points)*spectral_density_RR_2(j/disc_points)/disc_points;
    end
    variance_SGD_RR = variance_SGD_RR + variance_RR_i;
end
variance_SGD_RR
mean_MSE_SGD_RR_approx

for i=1:d
    evec = eigenvectors(:,i);
    sigma2 = evec.'*variance_A*evec
end

clear iterates_centered_SGD iterates_centered_SGD_approx 
clear iterates_centered_SGD_RR iterates_centered_SGD_RR_approx
clear iterates_centered_SGD_SO iterates_centered_SGD_SO_approx 
clear iterates_SGD iterates_SGD_approx 
clear iterates_SGD_RR iterates_SGD_RR_approx
clear iterates_SGD_SO iterates_SGD_SO_approx