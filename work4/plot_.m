clc; clear;

F = load('1.mat');
loadrun;
figure(1);
plotrun;
lambda1 = lambda;
e1 = e;

F = load('2.mat');
loadrun;
figure(2);
plotrun;
lambda2 = lambda;
e2 = e;

F = load('3.mat');
loadrun;
figure(3);
plotrun;
lambda3 = lambda;
e3 = e;

F = load('4.mat');
loadrun;
figure(4);
plotrun;
lambda4 = lambda;
e4 = e;

figure(5);
plot(time, lambda1, time, lambda2, time, lambda3, time, lambda4, time, lambda_opt);
xlabel('time [s]')
ylabel('travel [rad]')
legend('$\lambda_1$', '$\lambda_2$', '$\lambda_3$', '$\lambda_4$', '$\lambda_{opt}$', "Interpreter", "latex")

figure(6);
plot(time, e1, time, e2, time, e3, time, e4, time, e_opt);
xlabel('time [s]')
ylabel('elevation [rad]')
legend('$e_1$', '$e_2$', '$e_3$', '$e_4$', '$e_{opt}$', "Interpreter", "latex")
%}