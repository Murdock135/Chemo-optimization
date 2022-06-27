function dydt =ODEBryson(t,y,u)

global p r b a g s m myu ganma h KE KT UMAX UMIN Tmax W1 W2

u=u(floor(t+1));

dydt(1)= s- myu*y(1) + ((p*y(1)*y(2))/(h+y(2)))-m*y(1)*y(2)-KE*y(3)*y(1); % E(T)
dydt(2)=r*y(2)*(1-b*y(2))-(a*y(1)*y(2))/(y(2)+g)-KT*y(3)*y(2); %T(T)
dydt(3)= -ganma*y(3)+u; %M(T)
dydt(4)= W1*u+W2*y(2);
dydt=dydt';

end
