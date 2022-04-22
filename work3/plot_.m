clc; clear;

F = load('1.mat');
loadrun;
%figure(1);
%plotrun;
lambda1 = lambda;
p1 = p;

F = load('2.mat');
loadrun;
%figure(2);
%plotrun;
lambda2 = lambda;
p2 = p;

F = load('3.mat');
loadrun;
%figure(3);
%plotrun;
lambda3 = lambda;
p3 = p;

F = load('4.mat');
loadrun;
%figure(4);
%plotrun;
lambda4 = lambda;
p4 = p;

figure(5);
plot(time, lambda2, time, lambda3, time, lambda4, time, lambda1, time, lambda_opt);
xlabel('time [s]')
ylabel('travel [rad]')
legend('$\lambda_2$', '$\lambda_3$', '$\lambda_4$', '$\lambda_1$', '$\lambda_{opt}$', "Interpreter", "latex")

figure(6);
plot(time, p2, time, p3, time, p4);
xlabel('time [s]')
ylabel('pitch [rad]')
legend('$p_2$', '$p_3$', '$p_4$', "Interpreter", "latex")