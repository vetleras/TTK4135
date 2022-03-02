clc; clear;
T1 = load('1.mat').ans;
T2 = load('2.mat').ans;
T3 = load('3.mat').ans;

plot(T1(1,:), T1(2,:), T1(1,:), T2(2,:), T1(1,:), T3(2,:))
xlabel('time')
ylabel('travel')
legend('0.1', '1','10')