%Clearning the file 

clear all
clear figure
clc
format long

%% values of the parameters 

global p r b a g s m myu ganma h KE KT UMAX UMIN Tmax W1 W2
p = 0.015;
r = 4.31*(10^(-3)); %rate of tumor growth
b = 10^(-9); %capacity of the tumor cell
a = 3.41*(10^(-10)); % parameter of cancer cleanup 
g = 10^5; %half-saturation for cancer cleanup 
s = 1.2*(10^4); %growth rate of normal/effector cells
m = 2*(10^(-11));
myu = 4.12*(10^(-2));
ganma = 0.9;
h = 2.02*10;
KE = 6.00*(10^(-1));
KT = 8.00*(10^(-1));
UMAX = 10;% 1
UMIN = 0;
Tmax = 30;
W1=0.000001 ;
W2=1-W1;

%% initializations
options1=odeset('RelTol',1e-4,'AbsTol',1e-4); %analyze time frame
tmax=Tmax;
levels=3000;
dicretization=tmax/levels;
max_iter=30;%number of looping
tspan1=0:dicretization:tmax;%0~tmax in discretization level
len=length(tspan1)
u=10*ones(1,len); % initialize the input control 
count=1;
dhdu=zeros(1,len);
tspan2=tmax:-dicretization:0;%tmax~0 in -ve discretization level

%% this function to find time boundary for the singular arc 
c=s_arc(30000,40000,0,tspan1,tspan2,options1,u);

%% ODE equations 
for j=1:max_iter % 1 iteration
j;
[T,Y]=ode45(@(t,y) ODEBryson(t,y,u),tspan1,[30000 40000 0 0],options1);
%pause
y1t=[T'];
y1=Y(:,1);%y1=E(t) (Effector cells) 
y2=Y(:,2);%y2=T(t) (Tumour cells)
y3=Y(:,3);%y3=M(t) (drug concentration in body)
J=Y(end,4);

%% %costate equations 
%function dydt = CostateBryson(t,x,x1t,x1)
[T1,cox]=ode45(@(t,x)  CostateBryson(t,x,y1,y2,y3,y1t,u),tspan2,[0, 0, 0],options1);

cox1=flipud(cox(:,1));
cox2=flipud(cox(:,2));
cox3=flipud(cox(:,3));
T1=flipud(T1);

%% 
cox11=cox(:,1);
cox22=cox(:,2);
cox33=cox(:,3);

t1=c(1,1);
t2=c(2,1);
E = y1;
T = y2;
M = y3;
t = y1t;
u = length(y1);
%% 

%____________________________________________________________________
% Here is the place for calling singular arc. dhdu > 0 u=umin, dhdu<0 u=UMAX;
imax= length(cox3);

for i=1:imax
    dhdu(i) = (1*W1)+cox3(i);
    dhdu(i) = round(dhdu(i),2);
    
    if (t(i)<t1 && dhdu(i)<0)
        u(i) = UMAX/2.25;
        
    elseif(t(i)<t1 && dhdu(i)>0)
        u(i) = UMIN;
    else        
        if (t(i)>t2 && dhdu(i)>0)    
            u(i) = UMIN;
        elseif(t(i)>t2 && dhdu(i)==0)
            u(i) = 0.2188*t(i)^(0.0245);
        elseif(t(i)>t2 && dhdu(i)<0)
            u(i) = UMAX/2.25;
        else          
           u(i) = 1.7;  
        end
    end
end


%% Try state constraint function

[T,Y]=ode45(@(t,y) try_st(t,y,c,y1t,y1,y2,y3,u),tspan1,[30000, 40000, 0, 0],options1);

y1t=[T'];
yy1=Y(:,1);
yy2=Y(:,2);
yy3=Y(:,3);
yy4=Y(:,4);
yy4(end)

y11=flipud(yy1);
y22=flipud(yy2);
y33=flipud(yy3);
y1tt=flipud(y1t);

n = length(y1);

y1 = y11;
y2 = y22;
y3 = y33;

L1 = cox11;
L2 = cox22;
L3 = cox33;

%% The Equation of the langragian constant 

for i = imax:-1:1
    if (t(i)>t2) 
        n(i) = 0;
    else 
        if (t(i)<t1)   
            n(i) = 0;       
        else
           P1 = -(L1(i)*p*y2(i)/(h+y2(i)));
           P2 = (L1(i)*m*y2(i));
           P3 = (L1(i)*myu);
           P4 = (L1(i)*KE*y3(i));
           P5 = (L2(i)*a*y2(i)/(y2(i)+g));
           
           LL2 = (-1*W2)+((L1(i)*p*y1(i)*y2(i))/(h+y2(i))^2)-((L1(i)*p*y1(i))/(h+y2(i)))+(L1(i)*m*y1(i))-(L2(i)*r)+(2*L2(i)*r*b*y2(i))-((L2(i)*a*y1(i)*y2(i))/(y2(i)+g)^2)+((L2(i)*a*y1(i))/(y2(i)+g))+(L2(i)*KT*y3(i));
           LL3 = (L1(i)*KE*y1(i))+(L2(i)*KT*y2(i))+(L3(i)*ganma);
           
           EE = s+p*((y1(i)*y2(i))/(h+y2(i)))-(m*y1(i)*y2(i))-(myu*y1(i))-(KE*y3(i)*y1(i));
           TT = r*y2(i)*(1-b*y2(i))-a*((y1(i)*y2(i))/(y2(i)+g))-(KT*y3(i)*y2(i));
           MM = -ganma*y3(i)+u(i);
           
           P6 = -KT*(LL2*y2(i)+TT*L2(i));
           P7 = -ganma*LL3;
           P8 = (P6+P7)/KE;
           P9 = -EE*L1(i);
           P10 = (P8+P9)/y1(i);
           
           n(i) = P1+P2+P3+P4+P5-P10;
        end
    end
end


%% backward integration to find the costates, time period [0 30] 
%[T1,X1]=ode45(@(t,y) try_Costate_constraints(t,y,y11,y22,y33,y1tt,cox11,cox22,cox33,c,u),tspan2,[0,0,0],options);
[T1,X1]=ode45(@(t,y) try_co(t,y,y11,y22,y33,y1tt,n),tspan2,[0,0,0],options1);

cox1=flipud(X1(:,1));
cox2=flipud(X1(:,2));
cox3=flipud(X1(:,3));
coxT=flipud(T1);

%
end

%__________________________________________________________________
%plottings 
[T,Y]=ode45(@(t,y) ODEBryson(t,y,u),tspan1,[30000 40000 0 0],options1);
%pause
y1=Y(:,1);%y1=E(t) (Effector cells) 
y2=Y(:,2);%y2=T(t) (Tumour cells)
y3=Y(:,3);%y3=M(t) (drug concentration in body)
J=Y(end,4);
x=1:len;

%plotting E
figure(1)
plot(T,yy1,'-','Linewidth',3);
legend('y1- white B cells');
pause(0.1)

%plotting T
figure(2)
plot(T,yy2,'-','Linewidth',3);
legend('y2- tumor cells');
pause(0.1)

%plotting M
figure(3)
plot(T,yy3,'-','Linewidth',3);
legend('y3- medicine');
pause(0.1)

%plotting the costate 
[T1,cox]=ode45(@(t,x)  CostateBryson(t,x,y1,y2,y3,y1t,u),tspan2,[0, 0, 0],options1);

cox1=flipud(cox(:,1));
cox2=flipud(cox(:,2));
cox3=flipud(cox(:,3));
x=1:len;

%plotting lambda 1= DH/dE
figure(4)
plot(T,cox1,'-','Linewidth',3);
legend('cox(1)');
pause(0.1)

%plotting lambda 2= DH/dT
figure(5)
plot(T,cox2,'-','Linewidth',3);
legend('cox(2)');
pause(0.1)

%plotting lambda 3= DH/dM
figure(6)
plot(T,cox3,'-','Linewidth',3);
legend('cox(3)');
pause(0.1)
  
%plot u 
dhdu(1:100)'
figure(7)
plot(T,u,'-','Linewidth',3);
legend('u');
pause(0.1)

%plot state equation dhdu
figure(8)
plot(T,dhdu,'-','Linewidth',3);
legend('dhdu');
pause(0.1)

figure(9)
plot(y1t, n,'-','LineWidth',2.0);
legend('n');
title('Lagrange Multiplier VS Time');
