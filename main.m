clc; clear all; close all;

%% parameters
global s mu  p m r b a g h y KE KT umax umin w1 w2
s = 1.2e4;
mu  = 4.12e-2;
p = 0.015;
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
tspan = linspace(0,1000,100000); % subdiving t=[t0 tf] into N subintervals
%% Initial conditions
E0= 30e3;
T0 = 40e3;
M0 = 0;
x0 = [E0 T0 M0]; % storing the initial values
u = 0*tspan; % no chemo drug in the beginning
% u=3
%% Solving the optimal control problem



% max_iter = 500;
% for i=1:max_iter
    [Tx,X] = ode45(@(t,x) states(t,x,u), tspan, x0); % getting the state trajectory for current u
    E = X(:,1);
    T = X(:,2);
    M = X(:,3);

    % solving the costate equations
    pfinal = [0 0 0]; % free-end problem
    tspan_reversed = flip(tspan);
    [Tp,P] = ode45(@(t,p) costateEq(), tspan_reversed, pfinal)

% end



