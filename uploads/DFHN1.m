clear all
clc
h=.01;fn=100;in=0;
t=in:h:fn;
n=((fn-in)/h);

for I=1:20
a=2.5;
b=2.5;
eps=.08;

X(1:n)=0;
Y(1:n)=0;
X(1) = -0.870;
Y(1) = -0.212;

for i=1:n 
    X(i+1) = X(i) + h*(X(i) - (X(i)^3 /3)- Y(i) + .1*I)  ;
    Y(i+1) = Y(i) + h*(eps*(a-Y(i)+ b* tanh(10*X(i))));
end

u=@(x,y)    x-(x.^3/3)-y+.1*I;
v=@(x,y)    eps*(a-x+b*tanh(10*x));

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
%set(gca,'white')
x1=linspace(-3,3,100);
f=inline(strcat('((a+b*tanh(10*x))-(x-(x.^3/3)+0.1*',int2str(I),')'));
fx=fzero(f,[0,1]);


y1=x1-(x1.^3/3)+.1*I;

y2=(a+b*tanh(10*x1));

subplot(2,1,1), plot(t,Y)
xlabel('t');
ylabel('v');
title('V\circ=V-(V^3/3)-W+I         W\circ=\phi(a+V-bW)      ','FontSize',14)
text(3,1.25,'V','color','b')
subplot(2,1,2), plot(x,y,'g',x1,y1,x1,y2,'r',fx,(fx+.7)/.8,'ko')
title('phase portrait')
text(-2.75,1,'V nullcline','color','b','FontSize',14)
text(1.5,2,'W nullcline','color','r','FontSize',14)
text(0,-1,'trajectory','color','g','FontSize',14)
axis([-3 3 -2 3])
xlabel('V','FontSize',14);
ylabel('W','FontSize',14);
%F(I) = getframe(gca)

pause(.5)



     %h = gcf;
     %M(I) = getframe(h,[5 5 480 380]);
    %  hold off
end
        %  movie(M);
     %movie2avi(M,'fhn','fps',1)
%play and save movie
