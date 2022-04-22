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
plot(time, e, time, e_opt);
xlabel('time [s]')
ylabel('elevation [rad]')
legend('$e$', '$e_{opt}$', "Interpreter", "latex")

subplot(4,1,4);
plot(time, e_dot, time, e_dot_opt);
xlabel('time [s]')
ylabel('elevation [rad/s]')
legend('$\dot e$', '$\dot e_{opt}$', "Interpreter", "latex")