%%
n1 = 20; n2 = 30; n3 = 20;
theta1 = randn(1, n1); rho1 = ones(1, n1);
cplx_theta1 = exp(1i*theta1);
theta2 = randn(1, n2); rho2 = ones(1, n2);
cplx_theta2 = exp(1i*theta2);
theta3 = rand(1, n3) * 2 * pi; rho3 = ones(1, n2);
cplx_theta3 = exp(1i*theta3);

%%
figure
compass(cplx_theta1);
hold on
compass(mean(cplx_theta1))

figure
compass(cplx_theta2);
hold on
compass(mean(cplx_theta2))

figure
compass(cplx_theta3);
hold on
compass(mean(cplx_theta3))

figure
compass([mean(cplx_theta1),mean(cplx_theta2)])

%%
% resample n1 and n2 distributions
fullDist = [theta1, theta2];
testStat = abs(mean(cplx_theta1) - mean(cplx_theta2))

nBoot = 1000;
bootStat = zeros(1, nBoots);
for iBoot = 1 : nBoot
    bootPerm = randperm(length(fullDist));
    
    cplx_resamp1 = exp(1i * fullDist(bootPerm(1:n1)));
    cplx_resamp2 = exp(1i * fullDist(bootPerm(n1+1:end)));
    bootStat(iBoot) = abs(mean(cplx_resamp1) - mean(cplx_resamp2));
    
end

p = 1 - (find(bootStat > testStat, 1, 'first') / nBoot)
%%
figure
hist(bootStat)