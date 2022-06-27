function solns = state_eqns(t,x,u)
% x(1)=E
% x(2)=T
% x(3)=M
global s mu  k m r b a g h y KE KT umax umin w1 w2
t;
u=u(floor(t+1));
u;
dx = zeros(3,1);
dx(1) = s-mu*x(1)+k*((x(1)*x(2))/(h+x(2)))-m*x(1)*x(2)-KE*x(3)*x(1);
dx(2) = r*x(2)*(1-b*x(2))-a*((x(1)*x(2))/(x(2)*g))-KT*x(3)*x(2);
dx(3) = -y*x(3)+u;
 J = w1*u+w2*x(2)

%solns = dx;
solns = [dx;J];