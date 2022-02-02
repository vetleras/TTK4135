%%%%%%%%%%%%%%%%%%%%%LABDAG 4
% FOR HELICOPTER NR 1-2
% This file contains the initialization for the helicopter assignment in
% the course TTK4115. Run this file before you execute QuaRC_ -> Build 
% to build the file heli_q8.mdl.

% Oppdatert høsten 2006 av Jostein Bakkeheim
% Oppdatert høsten 2008 av Arnfinn Aas Eielsen
% Oppdatert høsten 2009 av Jonathan Ronen
% Updated fall 2010, Dominik Breu
% Updated fall 2013, Mark Haring
% Updated spring 2015, Mark Haring


%%%%%%%%%%% Calibration of the encoder and the hardware for the specific
%%%%%%%%%%% helicopter
Joystick_gain_x = 1;
Joystick_gain_y = 1;


%%%%%%%%%%% Physical constants
g = 9.81; % gravitational constant [m/s^2]
l_c = 0.40; % distance elevation axis to counterweight [m]
l_h = 0.66; % distance elevation axis to helicopter head [m]
l_p = 0.175; % distance pitch axis to motor [m]
m_c = 1.92; % Counterweight mass [kg]
m_p = 0.65; % Motor mass [kg]
V_s0 = 5.7;
K_f = - (g*((m_c*l_c)-2*m_p*l_h))/(V_s0*l_h);
K_1 = K_f / (2 * m_p * l_p);
%K_1 = 0.941;
K_2 = (K_f * l_h) / ((m_c * (l_c)^2) + (2*m_p*(l_h)^2));
% K_3 = -(K_f*l_h*g)*(l_c*m_c-2*l_h*m_p)/((m_c*((l_c)^2))+2*m_p(((l_h)^2)+((l_p)^2)));

K_3 = (g*(m_c*l_c-2*m_p*l_h))/(m_c*(l_c).^2+2*m_p*((l_h).^2+(l_p).^2));

w = 1;
Cetta = 1.5;
K_pp = w^2/K_1;
K_pd = 2*Cetta*w/K_1;

%%%%%%%%%%%%%% Labdag 1 %%%%%%%%%%%%%%%%%
K_pp = lambda^2 / K_1;
K_pd = sqrt(4*K_pp/K_1);
s=(4*K_1*K_pp)/(K_pd*K_pd);

%%%%%%%%%%%%%%%%%%  Labdag 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[0 1 0; 0 0 0 ; 0 0 0];
B=[0 0; 0 K_1; K_2 0];
C=[1 0 0; 0 0 1];

Q=[100 0 0; 0 30 0; 0 0 100];
R=[2 0; 0 2];


K=lqr(A,B,Q,R);
F=inv(C*(inv(B*K-A))*B);
%F=[k_11 k_13; k_21 k_23]
% de to siste er for integralaffekten, vokser over tid
%bare forholdet mellom dem som har noe å si
%Nå har vi realtiv lave integraleffekt (de to siste verdiene)Vi velger
%verdier etter hvilken av de den ønsker å prioritere av p p_dot og e_dot.
%P_dot har vi ikke - den er 0, men stor verdi på denne vil hillikopteret
%motsette seg endringer. 

%%%%%%%%%%%  Integral LQR %%%%%%%%%%%%%%%%%%%%%%%%%%%
A_I = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0;1 0 0 0 0; 0 0 1 0 0];
B_I = [ 0 0;0 K_1 ;K_2 0; 0 0; 0 0];
G_I =[0 0; 0 0; 0 0; -1 0; 0 -1];
Q_I = diag([50 30 90 2 2]);  

R_I = diag([5 5]);  % forstyrrelse /error
C_I=[1 0 0 0 0; 0 0 1 0 0];
K_I = lqr(A_I,B_I,Q_I,R_I);
F_I = [K_I(1,1) K_I(1,3); K_I(2,1) K_I(2,3)];

%I = diag([1 1 1 1 1]);
%I=C_I*(inv(B_I*K_I-A_I))*(B_I*F_I+G_I);
%F_I=inv(C_I*(inv(B_I*K_I-A_I))*B_I);

%bruker verdier fra integraleffet lab 2 til å finne K verdier og da også

%%%%%%%%%%%%%%%%%% Labdag 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_lab3 = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 1 0;0 0 0 0 0; K_3 0 0 0 0];
B_lab3 = [ 0 0;0 K_1 ;0 0; K_2 0; 0 0];
%C_lab3 =[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]; 
C_lab3 =[1 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1]; % tar bort p_dot og e_dot

%r_C =rank(C_lab3)
K_I = lqr(A_I,B_I,Q_I,R);
F_lab3 = [K_I(1,1) K_I(1,3); K_I(2,1) K_I(2,3)];

%poler til systemet
Poles_cl = eigs(A-B*K); % poles av det lukket systemet, lab 2 systemet


P = [-10,-10,-10,-10,-10];
L = place(A_lab3.',C_lab3.',P).';

%Poles_test = eigs(A_I-B_I*K_I);
%poler til lukket system, tommelfingeregel
%L = transpose(place(transpose(A_I),transpose(C_I),P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  lab 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vi laget en ny matrise for å måle A og estimere

Est = load('A.mat');
mal = Est.lab2_o2_01(2:7,1500:end);

% estimerte R_d verdier for systemet både når hellikopteret er stille
% og når hellikoppteret flyr i lineriseringspunktet
R_d = cov(mal')


R_d =   [  0.0027    0.0004    0.0003   -0.0023   -0.0010;
           0.0004    0.0278    0.0355   -0.0064    0.0032;
           0.0003    0.0355    0.0541   -0.0102    0.0048;
          -0.0023   -0.0064   -0.0102    0.0183   -0.0006;
          -0.0010    0.0032    0.0048   -0.0006    0.0081];


% R_d =   [0.0047   -0.0021    0.0001   -0.0031   -0.0289   -0.0042;
%        -0.0021    0.0238    0.0291   -0.0077    0.0166    0.0040;
%         0.0001    0.0291    0.0461   -0.0094    0.0074    0.0023;
%        -0.0031   -0.0077   -0.0094    0.0157   -0.0031    0.0017;
%        -0.0289    0.0166    0.0074   -0.0031    3.2798    0.2577;
%        -0.0042    0.0040    0.0023    0.0017    0.2577    0.0238];

% R_d =   [0.0047   -0.0021    0.0001   -0.0031   -0.0042;
%        -0.0021    0.0238    0.0291   -0.0077    0.0040;
%         0.0001    0.0291    0.0461   -0.0094   0.0023;
%        -0.0031   -0.0077   -0.0094    0.0157   0.0017;
%        -0.0042    0.0040    0.0023    0.0017  0.0238];
   
%R_d_stille = 1.0e-05 * [  0.1280   -0.0015    0.0004    0.0123         0    0.0015;
%                             -0.0015    0.0429   -0.0073    0.0001         0   -0.0083;
%                             0.0004   -0.0073    0.2663   -0.0118         0    0.0105;
%                             0.0123    0.0001   -0.0118    0.2307         0   -0.0167;
%                             0         0         0         0             0         0;
%                             0.0015   -0.0083    0.0105   -0.0167         0    0.0176];
       

% p,pd,e,ed,t,td
  %Q_d = diag([10 100 35 1 20 3])*10^(-4);%*10^(-4); %diag([0 0 0 0 0 0]);
  
  Q_d = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])^(-3);
  
%   Q_d= [0.0047   -0.0021    0.0001   -0.0031   -0.0289   -0.0042;
%        -0.0021    0.0238    0.0291   -0.0077    0.0166    0.0040;
%        0.0001    0.0291    0.0461   -0.0094    0.0074    0.0023;
%         -0.0031   -0.0077   -0.0094    0.0157   -0.0031    0.0017;
%         -0.0289    0.0166    0.0074   -0.0031    3.2798    0.2577;
%         -0.0042    0.0040    0.0023    0.0017    0.2577    0.0238];
%   
   %
   P_init = eye(6);% %diag([0.0122 0.0114 0.0234 0.0151 63.6004 0.01646]); %eye(6);%zeros(6,6); 
   %her forteller vi hvor gode målingene vår er,
   %normalfordleing/standardavvik
   %
   X_init =  [-0.01; 0; -0.54; 0; 0; 0]; %([1 1 1 1 1 1])'; %
   %denne forteller hva vi tror verdiene vår er. 
                        
   % nye matriser til systemet vi bruker nå alle tilstander
   % fordi t er ikke målbar med IMU
   
      
%    B_d = [ 0 0;
%            K_1 0;
%            0 0;
%            0 K_2;
%            0 0;
%            0 0];
% 
  
A_d = [ 0 1 0 0 0 0;
           0 0 0 0 0 0;
           0 0 0 1 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 1;
          K_3 0 0 0 0 0];
       
 B_d = [ 0 0;
         0 K_1;
         0 0;
         K_2 0;
         0 0;
         0 0];
       
     C_d = [ 1 0 0 0 0 0;
           0 1 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 1 0 0;
           0 0 0 0 0 1]; 
     
D_d = 0;
sys = ss(A_d, B_d, C_d, D_d); 
d_sys = c2d(sys, 0.002);
[A_d, B_d,C_d, D_d] = ssdata(d_sys);  
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot lab 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 T = load('test9.mat');  %messurments
 % test =  diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-3)
 %test 2 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-4) vi
 %kjørte litt roigere
 %test 3 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-10)
 %test 4 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-1)
% test 5 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(10); % 
%test 6 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(100); 
%test 7 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*0;
% test 8 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-3) med
% new data = 1
% test 9 = diag([0.0047 0.0238 0.0461 0.0157 3.2798 0.0238])*10^(-3) med
% new data = 0 alternerende  pichen er vedlig god
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









PORT = 4;

