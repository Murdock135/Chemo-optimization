function dp = costateEq(t,p,E,T,M,u,Tx)
% p provides the initial values for the costate variables, p(1),p(2) and
% p(3)
global s mu k m r b a g h y KE KT umax umin w1 w2


u = u(floor(t+1));
u; % watch u
t; % watch t
E = interp1(Tx,E,t); % Interpolate the data %set (y1t, y1) at times t
T = interp1(Tx,T,t); % Interpolate the data %set (y1t, y1) at times t
M = interp1(Tx,M,t); % Interpolate the data %set (y1t, y1) at times t
dp = zeros(3,1);
p;% watch p
dp(1) = p(1)*(mu-k*T/(h+T)+m*T+KE*M) + p(2)*a*T/(T+g);
dp(2) = -1+m*E-p(1)*k*E*h/(h+T)^2 + p(2)*(-r+r*b+a*E*g/(T+g)^2+KT*M);
dp(3) = p(1)*KE*E + p(2)*KT*T + p(3)*y;
