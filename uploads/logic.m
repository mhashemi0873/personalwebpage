clc
clear all
for A=[1.4,.3]
a=-1;b=50;h=.01;tau=1;y0=.1;
n=tau/h;
T=((b-a)/h);
x=a:h:b;
y(1:T)=0;
    y(1:n+1)=y0;
for i=n+1:T
    y(i+1)=y(i)+h*A*y(i)*(1-y(i-n));
end
plot(x,y)
hold on
end