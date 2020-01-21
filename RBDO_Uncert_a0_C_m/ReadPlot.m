test = dlmread('MC_15000Samples.txt');
C_samples = test(:, 1);
epsF_samples = test(:, 2);

n = 5000 - 1;
[fC, xiC] = ksdensity(C_samples(end-n:end));
[fepsR, xiepsR] = ksdensity(epsF_samples(end-n:end));

figure();
subplot(2,1,1)
plot(C_samples, 'k');
xlabel('Step')
ylabel('C')
title('Markov Chain of C')

subplot(2,1,2)
plot(epsF_samples, 'k', 'linewidth', 1.05);
xlabel('Step')
ylabel('\sigma_{\epsilon R}')
title('Markov Chain of \sigma_{\epsilon R}')