%% parameters

s = 1.2e4;
mu = 4.12e-2;
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
%% solving state equations

syms E(t) T(t) M(t) u(t)

eqnE = diff(E,t) == s-mu*E+p*((E*T)/(h+T))-m*E*T-KE*M*E;
condE = E(0) == 3e4;

eqnT = diff(T,t) == r*T*(1-b*T)-a*((E*T)/(T+g))-KT*M*T;
condT = T(0) == 4e4;

eqnM = diff(M,t) == -y*M+u;
condM = M(0) == 0;
condu = u(0) == 0;

eqns = [eqnE;eqnT;eqnM]



