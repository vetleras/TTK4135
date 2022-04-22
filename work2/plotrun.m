loadrun;
subplot(4,1,1);
plot(time, lambda, time, lambda_opt);
xlabel('time [s]')
ylabel('travel [rad]')
legend('$\lambda$', '$\lambda_{opt}$', "Interpreter", "latex")

subplot(4,1,2);
plot(time, lambda_dot, time, lambda_dot_opt);
xlabel('time [s]')
ylabel('travel [rad/s]')
legend('$\dot \lambda$', '$\dot \lambda_{opt}$', "Interpreter", "latex")

subplot(4,1,3);
plot(time, p, time, p_opt);
xlabel('time [s]')
ylabel('pitch [rad]')
legend('$p$', '$p_{opt}$', "Interpreter", "latex")

subplot(4,1,4);
plot(time, p_dot, time, p_dot_opt);
xlabel('time [s]')
ylabel('pitch [rad/s]')
legend('$\dot p$', '$\dot p_{opt}$', "Interpreter", "latex")