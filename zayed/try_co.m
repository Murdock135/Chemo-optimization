function  differentials=try_co(t,x,y11,y22,y33,y1tt,n)
    global p r b a g s m myu ganma h KE KT UMAX UMIN Tmax W1 W2
    
    y1 = interp1(y1tt,y11,t);
    y2 = interp1(y1tt,y22,t);
    y3 = interp1(y1tt,y33,t);
    
    cox1=x(1);
    cox2=x(2); 
    cox3=x(3);
    
    n = interp1(y1tt,n,t);
    
    %y1'- lamda1'
    differentials(1)=(cox1*myu)- (cox1*p*(y2/(h+y2))) + (cox1*m*y2)+(cox1*KE*y3)+(cox2*a*(y2/(y2+g)))-n;
    %y2'- lamda2'
    differentials(2)= -1*W2+(cox1*p*((y1*y2)/((h+y2)^2)))-(cox1*p*((y1)/(h+y2)))+(cox1*m*y1)-(cox2*r)+(2*cox2*r*b*y2)-(cox2*a*((y1*y2)/((g+y2)^2)))+(cox2*a*((y1)/(g+y2)))+(cox2*KT*y3);
    %y3'- lamda3'
    differentials(3)= (cox1*KE*y1)+(cox2*KT*y2)+(cox3*ganma);
    differentials=differentials';
end



