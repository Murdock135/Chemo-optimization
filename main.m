clc; clear all; close all;

%% parameters
global s mu  k m r b a g h y KE KT umax umin w1 w2
s = 1.2e4;
mu  = 4.12e-2;
k = 0.015;
m = 2e-11;
r = 4.31e-3;
b = 10^-9;
a = 3.41e-10;
g = 10^5;
h = 2.2;
y = 0.9;
KE = 6e-1;
KT = 8e-1;
umax = 1;
umin = 0;
w1 = 1;
w2 = 1;
%% Setting the window
interval = 0:5; % subdiving t=[t0 tf] into N subintervals
%% Initial conditions
E0= 30e3;
T0 = 40e3;
M0 = 0;
x0 = [E0 T0 M0]; % storing the initial values
u = 0*interval; % no chemo drug in the beginning
% u=3
%% Solving the optimal control problem



% max_iter = 500;
% for i=1:max_iter
    [Tx,X] = ode45(@(t,x) states(t,x,u), interval, [x0 0]); % getting the state trajectory for current u
    E = X(:,1);
    T = X(:,2);
    M = X(:,3);
    J1 = X(:,4);

    % solving the costate equations (backwards)
    pfinal = [0 0 0]; % free-end problem
    interval_reversed = flip(interval);

    [Tp,P_reverse] = ode45(@(t,p) costateEq(t,p,E,T,M,u,Tx), interval_reversed, pfinal);
    P = flip(P_reverse); % getting costate trajectory (forwards)
    p1 = P(:,1);
    p2 = P(:,2);
    p3 = P(:,3);
    % computing the cost at every time point/step
    u = u';
    J2 = u+T; 
    J = [J1 J2];
    % Computing dH/du
    H = u+T+p1.*(s-mu.*E+k.*((E.*T)./(h+T))-m.*E.*T-KE.*M.*E)...
        +p2.*(r.*T.*(1-b.*T)-a.*((E.*T)./(T.*g))-KT.*M.*T)...
        +p3.*(-y.*M+u);

% end

%% plot of state trajectory
figure(1)
plot(E)
hold on
plot(T)
hold on
plot(M)

legend('E','T','M')
title('state trajectory')

figure(2)
plot(J)
hold on
plot(H)
legend('J','dH/du')
title('cost, H and dH/du')




