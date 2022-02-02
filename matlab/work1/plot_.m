T = load('test9.mat');
 time = T.lab2_o2_01(1,:);
 %ref = T.lab2_o2_01(2,:);

%%%%%%%%%%%%%%%%%%%% encoder rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = T.lab2_o2_01(2,:);
t_dot = T.lab2_o2_01(3,:);
p = T.lab2_o2_01(4,:);
p_dot = T.lab2_o2_01(5,:);
e = T.lab2_o2_01(6,:);
e_dot = T.lab2_o2_01(7,:);

%%%%%%%%%%%%%%%%%%% estimator rates %%%%%%%%%%%%%%%%%%%%%%%%
p_hat = T.lab2_o2_01(8,:);
p_hat_dot = T.lab2_o2_01(9,:);
e_hat = T.lab2_o2_01(10,:);
e_hat_dot = T.lab2_o2_01(11,:);
t_hat = T.lab2_o2_01(12,:);
t_hat_dot = T.lab2_o2_01(13,:);


%%%%%%%%%%%%%%%%% euler rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 p_eu = T.lab2_o2_01(14,:);
 p_eu_dot = T.lab2_o2_01(15,:); %gt_p
 e_eu = T.lab2_o2_01(16,:);
 e_eu_dot = T.lab2_o2_01(17,:); %gt_e
 t_eu_dot = T.lab2_o2_01(18,:); %gt_t
 
%%%%%%%%% plot
% BLÅ = euler/y
% RØD = x_hat
% GUL encoder 

%legend
subplot(3,2,1)
plot(time, p_eu, time, p_hat,time, p)
xlabel('time')
ylabel('pitch')
title('pitch')
legend('x_{eu}', 'x_{hat}','x')

subplot(3,2,2)
plot( time, p_eu_dot,time,p_hat_dot,time, p_dot)
xlabel('time')
ylabel('pitch rate')
title('pitch rate')
legend('x_{eu}', 'x_{hat}','x')

subplot(3,2,3)
plot( time, e_eu, time, e_hat,time, e)
xlabel('time')
ylabel('elevation')
title('elevation')
legend('x_{eu}', 'x_{hat}','x')

subplot(3,2,4)
plot( time, e_eu_dot, time, e_hat_dot, time, e_dot)
xlabel('time')
ylabel('elevation rate')
title('elevation rate')
legend('x_{eu}', 'x_{hat}','x')

subplot(3,2,5)
plot(time, t, time, t_hat ,time, t) % time, t_eu
xlabel('time')
ylabel('travel')
title('travel')
legend('x_{eu}', 'x_{hat}','x')

subplot(3,2,6)
plot( time,- t_eu_dot, time, -t_hat_dot,time, t_dot) % er det feil fortegn?
xlabel('time')
ylabel('travel rate')
title('travel rate')
legend('x_{eu}', 'x_{hat}','x')

% 
% %%%%%%%%%%%%%%%%%%%%% joystick / controller %%%%%%%%%%%%%%%%%%%%%
% e_c_dot = T.lab2_o2_01(18,:);
% p_c = T.lab2_o2_01(19,:);
% 
% %%%%%%%%%%%%%%%%%%%%%% acc %%%%%%%%%%%%%%%%%%%%%%%%%
% a_x = T.lab2_o2_01(20,:);
% a_y = T.lab2_o2_01(21,:);
% a_z = T.lab2_o2_01(22,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%% gyro %%%%%%%%%%%%%%%%%%%%%%%%
% w_x = T.lab2_o2_01(23,:);
% w_y = T.lab2_o2_01(24,:);
% w_z = T.lab2_o2_01(25,:);


%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(time, p, time, p_eu+0.0464);
% 
% plot(time, p-p_eu);
% %R_v_p = mean(p-p_eu)   % gjennomsnittt av R_v (målefeil) på pitch 0.0463
% %p_cov = cov(p,p_eu)
% e_cov = cov(e, e_eu)
% %R_d = [];  % 5x5 matrise 
% 
% %plot(time, e-e_eu);
% o = mean(e-e_eu)
% 
% %plot(time, p_cov);
% %a = autocorr(p-p_eu)



%plot(time, p,time, p_eu,  time, p_hat);
%plot(time, p_dot, time, p_eu_dot, time, p_hat_dot);
%plot(time, e, time, e_eu, time, e_hat);
%plot(time, e_dot, time, e_eu_dot, time, e_hat_dot);
%plot(time, -t_dot, time, t_eu_dot, time, t_hat_dot);


%plot(time,elev, time, e_dot, time, e_c_dot);
%plot(time,pitch,time,elev,time,travel);
%plot(time,p_rate,time,e_rate,time,t_rate);


% plot(time,gt_t,time,-t_rate);


%plot(time, t_eu,time,t);
%plot(time, e_eu_dot,time,e_dot);
%plot(time, t_eu_dot,time,-t_dot);

%plot(time,w_x,time,p_dot);
%plot(time,w_y,time,e_dot);
%plot(time,w_z,time,t_dot);
%plot(



%plot(time,a_x,time,a_y,time,a_z);
%plot(time,w_x,time,w_y,time,-w_z);
%plot(time,a_z,time,w_z);
%ploplot(time,w_z,);
