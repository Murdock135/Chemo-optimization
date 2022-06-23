function time=s_arc(E0,T0,M0,tspan1,tspan2,options,u)

    global p r b a g s m myu ganma h KE KT UMAX UMIN Tmax W1 W2
    t1=30;
    t2=30;
    
    [T,Y]=ode15s(@(t,y) ODEB(t,y,u),tspan1,[E0 T0 M0],options);

    y1t=[T'];
    y1=Y(:,1);%y1=E(t) (Effector cells) 
    y2=Y(:,2);%y2=T(t) (Tumour cells)
    y3=Y(:,3);%y3=M(t) (drug concentration in body)
 
    max=length(y1);

    for i = 1:max
        if(y1(i)<10000)
            t1 = y1t(i-1);
            break
        end
    end
    
    [T,Y]=ode15s(@(t,y) ODEBryson2(t,y,u,y1t,y1,y2,y3),tspan1,[E0 T0 M0],options);

    y1t=flipud([T']);
    y1=flipud(Y(:,1));%y1=E(t) (Effector cells) 
    y2=flipud(Y(:,2));%y2=T(t) (Tumour cells)
    y3=flipud(Y(:,3));%y3=M(t) (drug concentration in body)
   
    [T1,cox]=ode15s(@(t,y) CostateBryson(t,y,y1,y2,y3,y1t,u),tspan2,[0, 0, 0],options);

    cox1=flipud(cox(:,1));
    cox2=flipud(cox(:,2));
    cox3=flipud(cox(:,3));
    T1=flipud(T1);
   
    for i=1:max
        v=W1+cox3(i);
        if(abs(v)<0.000001)   
            t2=T1(i);
        end
    end
    
    time=[t1;t2];

function dydt = ODEB(t,y,u)

    u=u(floor(t+1));
    dydt(1) = s+p*((y(1)*y(2))/(h+y(2)))-(m*y(1)*y(2))-(myu*y(1))-(KE*y(3)*y(1));
    dydt(2) = r*y(2)*(1-b*y(2))-a*((y(1)*y(2))/(y(2)+g))-(KT*y(3)*y(2));
    dydt(3) = -ganma*y(3)+u;
    dydt=dydt';
end

function dydt = ODEBryson2(t,y,u,y1t,y1,y2,y3)

    y1 = interp1(y1t,y1,t);
    y2 = interp1(y1t,y2,t);
    y3 = interp1(y1t,y3,t);
    u=u(floor(t+1));
   
    if (t>t1)
        u = 1.7;
       
    end
        
    dydt(1) = s+p*((y(1)*y(2))/(h+y(2)))-(m*y(1)*y(2))-(myu*y(1))-(KE*y(3)*y(1));
    dydt(2) = r*y(2)*(1-b*y(2))-a*((y(1)*y(2))/(y(2)+g))-(KT*y(3)*y(2));
    dydt(3) = -ganma*y(3)+u;
    dydt=dydt';

end

function differentials = CostateBryson(t,x,y1,y2,y3,y1t,u) %dH/dy = y'

    y1 = interp1(y1t,y1,t); 
    y2 = interp1(y1t,y2,t); 
    y3 = interp1(y1t,y3,t); 

    cox1=x(1);
    cox2=x(2); 
    cox3=x(3);

    %y1'- lamda1'
    differentials(1)=(cox1*myu)- (cox1*p*(y2/(h+y2))) + (cox1*m*y2)+(cox1*KE*y3)+(cox2*a*(y2/(y2+g)));
    %y2'- lamda2'
    differentials(2)= -1*W2+(cox1*p*((y1*y2)/((h+y2)^2)))-(cox1*p*((y1)/(h+y2)))+(cox1*m*y1)-(cox2*r)+(2*cox2*r*b*y2)-(cox2*a*((y1*y2)/((g+y2)^2)))+(cox2*a*((y1)/(g+y2)))+(cox2*KT*y3);
    %y3'- lamda3'
    differentials(3)= (cox1*KE*y1)+(cox2*KT*y2)+(cox3*ganma);
    differentials=differentials';
end

end