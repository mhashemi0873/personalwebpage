clear all
clc
h=.01;fn=20;in=0;
t=in:h:fn;
n=((fn-in)/h);

I=0;%input('input the value for current I=');
a=2.5;
b=2.5;
eps=.08;

X(1:n)=0;
Y(1:n)=0;
X(1) = 0;
Y(1) = 0;

for i=1:n 
    X(i+1) = X(i) + h*(X(i) - (X(i)^3 /3)- Y(i) + I)  ;
    Y(i+1) = Y(i) + h*(eps*(a -Y(i)+ b*tanh(10*X(i))));
end

u=@(x,y)    x-(x.^3/3)-y+I;
v=@(x,y)    eps*(a-y+b*tanh(10*x));

x0=1;y0=-1;

start=[x0,y0];
delt=.05;
nstep=10000;
x=zeros(1,nstep+1);
y=x;
x(1)=start(1);
y(1)=start(2);

for n=1:nstep
    x(n+1)=x(n)+delt*u(x(n),y(n));
    y(n+1)=y(n)+delt*v(x(n),y(n));
end;

x1=linspace(-3,3,100);

y1=x1-(x1.^3/3)+I;

y2=(a+b*tanh(10*x1));

subplot(2,2,1), plot(t,Y)
xlabel('time');
ylabel('V');
title('I=.0     ')
subplot(2,2,3), plot(x,y,'g',x1,y1,x1,y2,'r')
title('phase portrait')
text(-2.5,6,'V nullcline','color','b')
text(.5,3.5,'W nullcline','color','r')
text(0,-2,'trajectory','color','g')
axis([-3 3 -5 8])
xlabel('V');
ylabel('W');


h=.01;fn=200;in=0;
t=in:h:fn;
n=((fn-in)/h);

I=1;%input('input the value for current I=');
a=2.5;
b=2.5;
eps=.08;

X(1:n)=0;
Y(1:n)=0;
X(1) = 0;
Y(1) = 0;

for i=1:n 
    X(i+1) = X(i) + h*(X(i) - (X(i)^3 /3)- Y(i) + I)  ;
    Y(i+1) = Y(i) + h*(eps*(a -Y(i)+ b*tanh(10*X(i))));
end

u=@(x,y)    x-(x.^3/3)-y+I;
v=@(x,y)    eps*(a-y+b*tanh(10*x));

x0=1;y0=-1;

start=[x0,y0];
delt=.05;
nstep=10000;
x=zeros(1,nstep+1);
y=x;
x(1)=start(1);
y(1)=start(2);

for n=1:nstep
    x(n+1)=x(n)+delt*u(x(n),y(n));
    y(n+1)=y(n)+delt*v(x(n),y(n));
end;

x1=linspace(-3,3,100);

y1=x1-(x1.^3/3)+I;

y2=(a+b*tanh(10*x1));

subplot(2,2,2), plot(t,Y)
xlabel('time');
ylabel('V');
title(        'I=1     ')
subplot(2,2,4), plot(x,y,'g',x1,y1,x1,y2,'r')
title('phase portrait')
text(-2.5,6,'V nullcline','color','b')
text(.5,3.5,'W nullcline','color','r')
text(0,-2,'trajectory','color','g')
axis([-3 3 -5 8])
xlabel('V');
ylabel('W');