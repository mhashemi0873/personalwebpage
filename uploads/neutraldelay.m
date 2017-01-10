clc
clear all
a=-1;b=5;h=.01;tau=1;
n=tau/h;
x=a:h:b;
T=((b-a)/h);
    y(1:n+1)=a:h:0;
for i=n+1:T
    y(i+1)=y(i)-y(i+1-n)+y(i-n);
end
plot(x,y)
