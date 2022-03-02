clc; clear;

F = load('1.mat');
loadrun;
subplot(2, 3, 1);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_1 = lambda;


F = load('2.mat');
loadrun;
subplot(2, 3, 2);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_2 = lambda;

F = load('3.mat');
loadrun;
subplot(2, 3, 3);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_3 = lambda;

F = load('4.mat');
loadrun;
subplot(2, 3, 4);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_4 = lambda;


F = load('5.mat');
loadrun;
subplot(2, 3, 5);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_5 = lambda;

F = load('6.mat');
loadrun;
subplot(2, 3, 6);
plot(time, p, time, p_opt);
xlabel('time')
ylabel('pitch')
legend('p', 'p_opt')
lambda_6 = lambda;

subplot(2, 3, 1);
plot(time, lambda_1, time, lambda_2, time, lambda_3, time, lambda_4, time, lambda_5, time, lambda_6);
xlabel('time')
ylabel('travel')
legend('1', '2', '3', '4', '5', '6');
