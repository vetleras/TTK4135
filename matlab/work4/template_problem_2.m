% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
init03; % Change this to the init file corresponding to your helicopter

A_c = [0 1 0 0 0 0;
       0 0 -K_2 0 0 0;
       0 0 0 1 0 0;
       0 0 -K_1*K_pp -K_1*K_pd 0 0;
       0 0 0 0 0 1;
       0 0 0 0 -K_3*K_ep -K_3*K_ed];

B_c = [0 0;
       0 0;
       0 0;
       K_1*K_pp 0;
       0 0;
       0 K_3*K_ep];

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A1 = eye(6) + delta_t*A_c;
B1 = delta_t*B_c;

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x0 = [pi 0 0 0 0 0]';           % Initial values

% Time horizon and initialization
N  = 40;                                  % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = ones(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds
ul 	    = -30*pi/180;                   % Lower bound on control
uu 	    = 30*pi/180;                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6
P1 = diag([10, 10]);                                % Weight on input
Q = 2.*gen_q(Q1,P1,N,M);                                  % Generate Q, hint: gen_q
f = @(x) 0.5*x'*Q*x;

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(N*mx, 1);             % Generate b
beq(1:6) = A1*x0;
opt = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',400000);
%% Solve QP problem with linear model
tic
[z,lambda] = fmincon(f, z0, [], [], Aeq, beq, vlb, vub, @(z)cons(z),opt); % hint: quadprog. Type 'doc quadprog' for more info 
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
%% Do more
u1  = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
u2  = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

t = 0:delta_t:delta_t*(length(u1)-1);

Q = eye(mx); R = diag([1, 2]); %test 1
Q = diag([1, 1, 1, 1, 10, 1]); R = diag([1, 1]); %test 2
Q = diag([1, 1, 1, 1, 1000, 1]); R = diag([1, 1]); %test 3
Q = diag([1, 1, 1, 1, 1000, 1]); R = diag([1, 10]); %test 4

Kt = dlqr(A1, B1, Q, R)';

tu = [t' u1 u2];
tx = [t' x1 x2 x3 x4 x5 x6];

%% Plotting


figure(2)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','mo'),grid
ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5,'mo'),grid
ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')
%}
