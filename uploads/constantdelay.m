clc
clear all
a=-1;b=25;h=.01;tau=1;
T=(b-a)/h;
n=tau/h;
x=a:h:b;
 y(1:T)=0;
 for i=1:n+1
     y(i)=1;
 end
for i=n+1:T
    y(i+1)=y(i)-h*y(i-n);
end
plot(x,y)
