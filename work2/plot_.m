clc; clear;

F = load('1.mat');
figure(1);
plotrun;
lambda1 = lambda;

F = load('2.mat');
figure(2);
plotrun;
lambda2 = lambda;

F = load('3.mat');
figure(3);
plotrun;
lambda3 = lambda;

figure(4);
plot(time, lambda1, time, lambda2, time, lambda3, time, lambda_opt);
xlabel('time [s]')
ylabel('travel [rad]')
legend('$\lambda_{q=0.1}$', '$\lambda_{q=1}$', '$\lambda_{q=10}$', '$\lambda_{opt q=10}$', "Interpreter", "latex")