clear all
clc
h=.01;fn=300;in=0;
t=in:h:fn;
n=((fn-in)/h);

I=.32;%input('input the value for current I=');
a=.7;
b=.8;
eps=.08;

X(1:n)=0;
Y(1:n)=0;
X(1) = -0.870;
Y(1) = -0.212;

for i=1:n 
    X(i+1) = X(i) + h*(X(i) - (X(i)^3 /3)- Y(i) + I)  ;
    Y(i+1) = Y(i) + h*(eps*(X(i) + a - b*Y(i)));
end

u=@(x,y)    x-(x.^3/3)-y+I;
v=@(x,y)    eps*(a+x-b*y);

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

y2=(x1+.7)/.8;

subplot(2,2,1), plot(t,Y)
xlabel('time');
ylabel('V');
title('I=.32      ')
%text(3,1.25,'V','color','b')
subplot(2,2,3), plot(x,y,'g',x1,y1,x1,y2,'r')
title('phase portrait')
text(-2.5,6,'V nullcline','color','b')
text(.5,4.5,'W nullcline','color','r')
%text(-2.5,6,'w=v-(v.^3/3)+I','color','b')
%text(2.,6,'w=(v+.7)/.8','color','r')
text(0,-1,'trajectory','color','g')
axis([-3 3 -3 10])
xlabel('V');
ylabel('W');


h=.01;fn=300;in=0;
t=in:h:fn;
n=((fn-in)/h);

I=.33;%input('input the value for current I=');
a=.7;
b=.8;
eps=.08;

X(1:n)=0;
Y(1:n)=0;
X(1) = -0.870;
Y(1) = -0.212;

for i=1:n 
    X(i+1) = X(i) + h*(X(i) - (X(i)^3 /3)- Y(i) + I)  ;
    Y(i+1) = Y(i) + h*(eps*(X(i) + a - b*Y(i)));
end

u=@(x,y)    x-(x.^3/3)-y+I;
v=@(x,y)    eps*(a+x-b*y);

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

y2=(x1+.7)/.8;

subplot(2,2,2), plot(t,Y)
xlabel('time');
ylabel('V');
title(        'I=.33      ')
%text(3,1.25,'V','color','b')
subplot(2,2,4), plot(x,y,'g',x1,y1,x1,y2,'r')
title('phase portrait')
%text(-2.5,6,'w=v-(v.^3/3)+I','color','b')
%text(2.,6,'w=(v+.7)/.8','color','r')
text(-2.5,6,'V nullcline','color','b')
text(.5,4.5,'W nullcline','color','r')
text(0,-1,'trajectory','color','g')
axis([-3 3 -3 10])
xlabel('V');
ylabel('W');